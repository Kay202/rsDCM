// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 Godfred Arhin and Nilotpal Sanyal

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

static inline double mom_log_mix(const double beta, const double tau0,
                                 const double tau1, const double pi) {
  const double log_f0 = R::dnorm(beta, 0.0, tau0, true);
  const double log_f1 = (beta == 0.0) ? R_NegInf :
    -0.5 * std::log(2.0 * std::acos(-1.0)) - 3.0 * std::log(tau1) +
    2.0 * std::log(std::fabs(beta)) - 0.5 * beta * beta / (tau1 * tau1);
  const double a = std::log1p(-pi) + log_f0;
  const double b = std::log(pi) + log_f1;
  if (a >= b) return a + std::log1p(std::exp(b - a));
  return b + std::log1p(std::exp(a - b));
}

static inline double mom_resp(const double beta, const double tau0,
                              const double tau1, const double pi) {
  const double log_f0 = R::dnorm(beta, 0.0, tau0, true);
  const double log_f1 = (beta == 0.0) ? R_NegInf :
    -0.5 * std::log(2.0 * std::acos(-1.0)) - 3.0 * std::log(tau1) +
    2.0 * std::log(std::fabs(beta)) - 0.5 * beta * beta / (tau1 * tau1);
  const double d = (std::log(pi) + log_f1) - (std::log1p(-pi) + log_f0);
  if (d >= 0.0) return 1.0 / (1.0 + std::exp(-d));
  const double ed = std::exp(d);
  return ed / (1.0 + ed);
}

static inline double t_log_density_residual(const double residual,
                                             const double inv_nu_sigma2,
                                             const double log_constant,
                                             const double exponent) {
  const double scaled_square = residual * residual * inv_nu_sigma2;
  return log_constant - exponent * std::log1p(scaled_square);
}

// [[Rcpp::export]]
List mom_parallel_tempering_cpp(NumericVector beta_init,
                                NumericVector y_star,
                                NumericMatrix X_star,
                                NumericVector inverse_temperatures,
                                NumericVector coordinate_sd,
                                double nu, double sigma2,
                                double tau0, double tau1, double pi,
                                int n_iter, int burn, int thin,
                                int n_chains, int n_temperatures,
                                double large_jump_prob, int swap_every) {
  const int q = beta_init.size();
  const int n = y_star.size();
  const int nt = n_temperatures;
  const int nkeep = (n_iter - burn) / thin;
  const double inv_nu_sigma2 = 1.0 / (nu * sigma2);
  const double log_constant = std::lgamma(0.5 * (nu + 1.0)) -
    std::lgamma(0.5 * nu) - 0.5 * std::log(nu * std::acos(-1.0) * sigma2);
  const double t_exponent = 0.5 * (nu + 1.0);
  Function sample_int = Environment::base_env()["sample.int"];

  List pip_draws(n_chains);
  NumericMatrix beta_means(n_chains, q);
  NumericMatrix local_acceptance(n_chains, nt);
  NumericMatrix large_jump_acceptance(n_chains, nt);
  NumericMatrix swap_acceptance(n_chains, nt - 1);
  NumericMatrix proposal_scales(n_chains, nt);

  for (int chain = 0; chain < n_chains; ++chain) {
    NumericMatrix states(nt, q);
    for (int t = 0; t < nt; ++t)
      for (int j = 0; j < q; ++j) states(t, j) = beta_init[j];

    for (int t = 0; t < nt; ++t) {
      if (chain > 0 || t > 0) {
        const double slab_probability = std::max(
          static_cast<double>(chain) / n_chains,
          static_cast<double>(t) / nt
        );
        std::vector<int> slab(q, 0);
        for (int j = 0; j < q; ++j)
          slab[j] = R::runif(0.0, 1.0) < slab_probability;
        for (int j = 0; j < q; ++j)
          if (slab[j]) states(t, j) += R::rnorm(0.0, tau1);
        for (int j = 0; j < q; ++j)
          if (!slab[j]) states(t, j) += R::rnorm(0.0, tau0);
      }
    }

    NumericMatrix residuals(n, nt);
    NumericMatrix prior_terms(nt, q);
    NumericVector log_likelihoods(nt), prior_sums(nt), log_posteriors(nt);
    for (int t = 0; t < nt; ++t) {
      for (int i = 0; i < n; ++i) {
        double fitted = 0.0;
        for (int j = 0; j < q; ++j) fitted += X_star(i, j) * states(t, j);
        residuals(i, t) = y_star[i] - fitted;
      }
      long double ll_acc = 0.0L;
      for (int i = 0; i < n; ++i) {
        ll_acc += static_cast<long double>(t_log_density_residual(
          residuals(i, t), inv_nu_sigma2, log_constant, t_exponent));
      }
      log_likelihoods[t] = static_cast<double>(ll_acc);
      double psum = 0.0;
      for (int j = 0; j < q; ++j) {
        prior_terms(t, j) = mom_log_mix(states(t, j), tau0, tau1, pi);
        psum += prior_terms(t, j);
      }
      prior_sums[t] = psum;
      log_posteriors[t] = log_likelihoods[t] + psum;
    }

    NumericVector proposal_scale(nt);
    std::fill(proposal_scale.begin(), proposal_scale.end(), 2.4);
    std::vector<int> accepted_postburn(nt, 0), proposed_postburn(nt, 0);
    std::vector<int> accepted_window(nt, 0), proposed_window(nt, 0);
    std::vector<int> jump_accepted(nt, 0), jump_proposed(nt, 0);
    std::vector<int> swap_accepted(nt - 1, 0), swap_proposed(nt - 1, 0);
    int swap_phase = 0;
    int kept = 0;
    NumericMatrix pip_chain(nkeep, q);
    NumericVector beta_sum(q);

    for (int iteration = 1; iteration <= n_iter; ++iteration) {
      for (int t = 0; t < nt; ++t) {
        IntegerVector order = sample_int(q);
        for (int pos = 0; pos < q; ++pos) {
          const int j = order[pos] - 1;
          const bool is_large_jump = R::runif(0.0, 1.0) < large_jump_prob;
          const double temp_scale = 1.0 / std::sqrt(inverse_temperatures[t]);
          const double step_sd = is_large_jump ? tau1 * temp_scale :
            proposal_scale[t] * coordinate_sd[j] * temp_scale;
          const double proposed_value = states(t, j) + R::rnorm(0.0, step_sd);
          const double delta = proposed_value - states(t, j);
          long double proposed_ll_acc = 0.0L;
          for (int i = 0; i < n; ++i) {
            const double proposed_resid = residuals(i, t) - X_star(i, j) * delta;
            proposed_ll_acc += static_cast<long double>(t_log_density_residual(
              proposed_resid, inv_nu_sigma2, log_constant, t_exponent));
          }
          const double proposed_log_likelihood = static_cast<double>(proposed_ll_acc);
          const double proposed_log_prior_j = mom_log_mix(proposed_value, tau0, tau1, pi);
          const double proposed_log_prior_sum = prior_sums[t] - prior_terms(t, j) + proposed_log_prior_j;
          const double proposed_log_posterior = proposed_log_likelihood + proposed_log_prior_sum;
          const double log_ratio = inverse_temperatures[t] *
            (proposed_log_posterior - log_posteriors[t]);
          const bool accepted = std::isfinite(proposed_log_posterior) &&
            std::log(R::runif(0.0, 1.0)) < log_ratio;
          if (accepted) {
            states(t, j) = proposed_value;
            for (int i = 0; i < n; ++i)
              residuals(i, t) -= X_star(i, j) * delta;
            log_likelihoods[t] = proposed_log_likelihood;
            prior_terms(t, j) = proposed_log_prior_j;
            prior_sums[t] = proposed_log_prior_sum;
            log_posteriors[t] = proposed_log_posterior;
          }

          if (iteration > burn) {
            ++proposed_postburn[t];
            if (accepted) ++accepted_postburn[t];
            if (is_large_jump) {
              ++jump_proposed[t];
              if (accepted) ++jump_accepted[t];
            }
          } else if (!is_large_jump) {
            ++proposed_window[t];
            if (accepted) ++accepted_window[t];
          }
        }
      }

      if (iteration <= burn && iteration % 50 == 0) {
        for (int t = 0; t < nt; ++t) {
          const double rate = static_cast<double>(accepted_window[t]) /
            std::max(proposed_window[t], 1);
          const double adj = std::max(-0.5, std::min(0.5, rate - 0.44));
          proposal_scale[t] *= std::exp(adj);
          accepted_window[t] = proposed_window[t] = 0;
        }
      }

      if (iteration % swap_every == 0) {
        swap_phase = 1 - swap_phase;
        const int first_pair = swap_phase == 0 ? 0 : 1;
        for (int lower = first_pair; lower < nt - 1; lower += 2) {
          const int upper = lower + 1;
          const double log_swap_ratio =
            (inverse_temperatures[lower] - inverse_temperatures[upper]) *
            (log_posteriors[upper] - log_posteriors[lower]);
          if (iteration > burn) ++swap_proposed[lower];
          if (std::log(R::runif(0.0, 1.0)) < log_swap_ratio) {
            for (int j = 0; j < q; ++j) std::swap(states(lower, j), states(upper, j));
            for (int i = 0; i < n; ++i) std::swap(residuals(i, lower), residuals(i, upper));
            for (int j = 0; j < q; ++j) std::swap(prior_terms(lower, j), prior_terms(upper, j));
            std::swap(log_likelihoods[lower], log_likelihoods[upper]);
            std::swap(prior_sums[lower], prior_sums[upper]);
            std::swap(log_posteriors[lower], log_posteriors[upper]);
            if (iteration > burn) ++swap_accepted[lower];
          }
        }
      }

      if (iteration > burn && (iteration - burn) % thin == 0) {
        for (int j = 0; j < q; ++j) {
          pip_chain(kept, j) = mom_resp(states(0, j), tau0, tau1, pi);
          beta_sum[j] += states(0, j);
        }
        ++kept;
      }
    }

    pip_draws[chain] = pip_chain;
    for (int j = 0; j < q; ++j) beta_means(chain, j) = beta_sum[j] / nkeep;
    for (int t = 0; t < nt; ++t) {
      local_acceptance(chain, t) = proposed_postburn[t] > 0 ?
        static_cast<double>(accepted_postburn[t]) / proposed_postburn[t] : NA_REAL;
      large_jump_acceptance(chain, t) = jump_proposed[t] > 0 ?
        static_cast<double>(jump_accepted[t]) / jump_proposed[t] : NA_REAL;
      proposal_scales(chain, t) = proposal_scale[t];
    }
    for (int t = 0; t < nt - 1; ++t)
      swap_acceptance(chain, t) = swap_proposed[t] > 0 ?
        static_cast<double>(swap_accepted[t]) / swap_proposed[t] : NA_REAL;
  }

  return List::create(
    _["pip_draws"] = pip_draws,
    _["beta_means"] = beta_means,
    _["local_acceptance"] = local_acceptance,
    _["large_jump_acceptance"] = large_jump_acceptance,
    _["swap_acceptance"] = swap_acceptance,
    _["proposal_scales"] = proposal_scales
  );
}

