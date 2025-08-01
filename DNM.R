# Define input parameters
observed_mutations <- NA  # Observed number of de novo mutations
total_population <- NA    # Total number of individuals in the study
mutation_rate <- NA       # Per gene, type-specific mutation rate (e.g., from Nguyen 2017)

# Calculate expected mutation rate per individual (times 2 for autosomal genes)
expected_mutation_rate <- 2 * mutation_rate

# Calculate the expected number of mutations
expected_mutations <- total_population * expected_mutation_rate

# Perform one-sided Poisson test
# ppois gives the probability of observing up to a certain number of events.
# For one-tailed test (greater side), subtract this from 1.
p_value <- 1 - ppois(observed_mutations - 1, expected_mutations)
pval <- ppois(observed_mutations - 1, lambda = expected_mutations, lower.tail = FALSE)

# Print the p-values (both should be identical)
print(p_value)
print(pval)


