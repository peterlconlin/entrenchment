# Entrenchment
## Supplementary Material for "Division of labor promotes the entrenchment of multicellularity"

**Authors:** Peter L. Conlin, Heather J. Goldsby, Eric Libby, Katherine G. Skocelas, William C. Ratcliff, Charles Ofria, and Benjamin Kerr

---

### Materials and Methods

#### Avida

Digital evolution experiments leverage the power of computers to ensure rapid generations (millions within days), perfect data collection (evolutionary lineages stored at fine temporal resolution), and high levels of experimental control (selective pressures toggled with configuration files). Here we use the software platform Avida [1] which has previously been employed to study the evolutionary origin of complex features [2], division of labor [3, 4], and adaptive radiation [5]. In this study, we focus on the entrenchment of higher-level units formed via fraternal transitions, where a higher-level unit (a digital "multicell") is composed of a set of related lower-level units (digital "cells"). We maintain a population of 1,000 organisms, where each organism can be either unicellular or multicellular. Each cell consists of a program (i.e., its genome), where the instructions encode metabolism, development and reproduction.

##### Time

The standard unit of time in Avida is an "update." On average, cells receive 30 CPU cycles every update. The other unit of time in this study is a "generation." Along a line of descent, every time an offspring is produced the generation value increments by one. Within an evolving population, generation is calculated at the level of cells in the following way. Every cell gets the generation value of its parent cell incremented by one. Time in terms of generations is simply the mean generation of all cells within the population at any instant. Note that one unit of generation time within a population will track the average population turnover time of unicellular organisms, but will only be a fraction of the turnover time for multicellular organisms.

##### Instructions and functions

The Avidian genome is composed of instructions that enable cells to acquire resources, communicate with other cells in a multicellular context, reproduce within the body, and determine their eligibility to produce a propagule to initiate a new organism. A complete list of the instructions used in our evolution experiments, separated by category, are provided in Tables 1-5.

Table 1: **No-operation instructions.** Instructions used in this study that have no direct effects when executed.

Table 2: **Hardware control instructions.** Instructions used in this study that control the reading/writing of the Avidian genome and handle memory.

Table 3: **Math instructions.** Instructions used in this study that execute mathematical operations.

Table 4: **Interaction instructions.** Instructions used in this study that enable communication and interaction among cells within a multicellular organism.

Table 5: **Biological instructions.** Instructions used in this study that underpin cellular differentiation, multicellular development, and reproduction.

The metabolism of digital cells involves the completion of logic functions, where the execution of a series of genomic instructions enables the cell to export the result of a logic operation on input bitstrings. The nine possible logic functions that can be performed are NOT, NAND, AND, ORNOT, OR, ANDNOT, NOR, XOR, and EQUALS. The execution of each

### Workload

<!--The dirty work paper has a really good explanation of workload. it's here: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001858#s3 in the supplement-->
A cell's workload is defined as:
<!--Cell workload can't be defined in the same way as it was done for the dirty work paper because here we have NAND as non-mutagenic (in addition to NOT). Perhaps FML_AND goes in the denominator now? I'm really not certain.-->
w = ∑(i ∈ F) Ni * (FMLi / FML_AND)
<!-- Modified from PLOS Biology paper. -->
<!-- Need to check that this is correct!-->

where F ≡ {NOT, NAND, AND, ORNOT, OR, ANDNOT, NOR, XOR, EQUALS}, Ni is the number of times function i is performed, FMLi is the mutagenic level of function i (as specified in Table~t:FML), and FML_AND is the base FML (i.e., the mutagenic level of the function AND).
We define propagule workload difference as the mean workload of propagule-ineligible cells minus the mean workload of propagule-eligible cells.

### Mutations

Mutations to the genome, in which one instruction is replaced by another randomly selected instruction, occur during two different processes: organismal reproduction and the execution of "dirty work" (logic tasks with FML>0).
Mutations during reproduction only occur during propagule production (not tissue accretion) with probability 0.01 per site. On average this results in one new mutation per propagule produced. Dirty work mutations occur immediately after a cell performs a task associated with a strictly positive function mutagen level. Each site in the genome of the performing cell experiences mutations with probability given by the FML of the task executed (Table~t:FML). More details behind the dirty work framework are provided in Goldsby et al. 2014~[goldsby2014evolutionary].

### Baseline evolution experiment

The ancestral organism is unicellular with a genome of 100 instructions. Execution of the ancestral genome initially enables only a few simple behaviors. Specifically, the ancestor completes the task NOT a single time, subsequently rotates counterclockwise and then tries to produce a propagule cell. This reproductive process fails since it only has collected ~5 units of resource. It repeats this execution loop until it has enough resources to produce a propagule (taking roughly 1800 updates for the first generation). The population has a maximum size of 1000 organisms and evolves for a total of 1 million updates with a reproductive mutation probability of 0.01 and a dirty work mutation probability defined by whatever mutagenic functions evolve in the organisms.

## Evolution of multicellularity

Within our system, multicellularity is generally beneficial to organisms as a result of the increased ability of multiple cells to rapidly access resources that can be devoted to propagule production. Here we explore how the efficiency of resource consumption affects the evolution of multicellularity. Within our original experiments, each cell was able to consume 5% of the available resources each time a task was performed regardless of the number of cells in the organism. Here, we observe what happens when we increase the efficiency of resource uptake from 5% to 8%. This efficiency is applied at the cell level and thus is available to cells within either a unicellular or multicellular context. As can be seen in Figure 1, increasing efficiency increases the percentage of replicates that evolve to be multicellular.

![Figure 1: Percent of evolutionary runs that evolve multicellularity as a function of resource acquisition efficiency.](figures/Figure_S1_Percent_multicellular_by_efficiency_level.pdf)

### Cost of multicellularity

We added an extrinsic cost of multicellularity for the purposes of testing how readily populations would revert to back to unicells. This cost was implemented as a time delay for producing a propagule from a multicellular organism. During this delay the organism and all its constituent cells are inert -- they cannot perform tasks, acquire resources, reproduce, or communicate. Additionally, the organism remains at risk for being replaced by a reproducing competitor.

To confirm the costly nature of this time delay, we reran the evolution of our initially unicellular digital populations varying this delay from 0 updates (as in our original treatment) to 512 updates. We varied this multicellularity cost from 0 updates (as in our original evolution treatment) to 512 updates, incrementing by powers of 2 (e.g., 0, 1, 2, 4, 8, 16, etc.). As expected, fewer replicates evolved multicells as this cost was increased (Figure 2).

![Figure 2: Percent of evolutionary runs that evolve multicellularity as a function of time delay cost.](figures/Percent_multi_by_cost.pdf)

### Entrenchment

In studies of molecular evolution, a focal substitution is said to be entrenched by subsequent substitutions if it becomes relatively more deleterious to revert as a result of those subsequent substitutions. Here, because we are examining the entrenchment of a *phenotype* (multicellularity), we require an entrenchment metric that integrates across all possible mutations that can cause reversion to unicellularity. Specifically, we assess the degree to which multicellularity has become entrenched, by performing the following steps. First, we create a new population with 1000 copies of the organism when it was born. We then add an externally imposed cost of multicellularity.  We run the population with this new multicellularity cost and observe if it reverts to unicellularity within 100 generations (assessed by measuring whether mean organism size has dropped below 2 cells after 100 generations of growth). We repeat this measurement with the following multicellularity cost values: 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048. To account for the stochastic nature of this assay, especially given the effect of dirty work on the organisms, we repeat the measurement three times for each organism, cost, and time point. We then compute a stability value for each genotype as the multicellular cost at which the probability of reversion to unicellularity $= 0.5$ using a logistic fit to the data (Figure 3). Finally, entrenchment is computed as the difference in stability value at the final and transition time points.

