functions {
    // Function to compute simple repression fold-change
    real log_fold_change(
        real logR,  // Log repressor copy number
        real eRA  // repressor binding energy
    ) {
        // Define fold-change
        real logfc = log(1 / (1 + exp(logR) / 4.6E6 * exp(-eRA)));

        return logfc;
    }
}

data{
    // Fold-change data
    real log_fc_exp[6, 4];  // log fold-change. col=operator, row=rep/cell
}

parameters {
    real<lower=1> logR[6];  // log repressors / cell
    // real<lower=0> sigma_logR[6];  // variance for log repressor distribution
    real<upper=0> eRA[4];  // repressor binding energies
    // real<lower=0> sigma_eRA[4]; // variance for binding energy
    real<lower=0> sigma_fc; // variance for binding energy
}

model {
    // Set Gaussian priors
    for (i in 1:6) {
        logR[i] ~ normal(1, 8.5);
    }
    for (i in 1:4) {
        eRA[i] ~ normal(0, 20);
    }
    sigma_fc ~ normal(0, 5);
    // Run inference model on data
    for (logr in 1:6) {
        for (era in 1:4) {
            // Compute theoretical fold-change
            real log_fc_thry = log_fold_change(logR[logr], eRA[era]);
            // Set Gaussian likelihood
            log_fc_exp[logr, era] ~ normal(log_fc_thry, sigma_fc);
        }
    }
}

// generated quantities {
//     // Generate posterior predictive samples for the fold-change
//     real log_fc_ppc[6, 4] = 
// }