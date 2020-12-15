# Seed setting for reproducible analyses

---

In DIANE, some of the available methods such as co-expression Clustering and Gene Regulatory Network Inference include random processes. This randomness is inherent to the chosen algorithms, and may produce different results for two identical runs.

To deal with this reality, that could get in the way of repeatable workflows, the seed, i.e the initial state of the random generator is set right before every analysis involving randomness so that two runs with identical settings will give the same outputs.

This seed is set here to a random value, but if you want to explore other possible outputs offered by randomness in DIANE, feel free to chose another seed by changing its value. You can also store the seed if you want to run the exact same analyses after exiting and reloading the app.
