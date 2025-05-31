# Create a numeric vector and store the vector as a variable called 'glengths'
glengths <- c(4.6, 3000, 50000)
glengths
# Create a character vector and store the vector as a variable called 'species'
species <- c("e.coli", "human", "corn")
species

combined <-  c(glengths, species)
combined

# Create a character vector and store the vector as a variable called 'expression'
expression <- c("low", "high", "medium", "high", "low", "medium", "high")

#Turn 'expression' vector into a factor
expression <- factor(expression)

#Exercises creating samplegroup with nine elements
samplegroup <- c("CTL", "KO", "OE", "CTL", "CTL", "OE", "KO", "KO", "OE")
samplegroup <- factor(samplegroup)

# Create a matrix with numbers from 1 to 9
my_matrix <- matrix(1:9, nrow = 3, ncol = 3)
my_matrix

# Create a data frame and store it as a variable called 'df'
df <- data.frame(species, glengths)
df

# Exercise
titles <- c("Catch-22", "Pride and Prejudice", "Nineteen Eighty Four")
pages <- c(453, 432,328)
favorite_books <- data.frame(titles, pages)
favorite_books


list1 <- list(species, df, number)
list1

list2 <- list(species, glengths, number)
list2

glengths <- c(glengths, 90)
glengths
glengths <- c(30, glengths)
glengths

sqrt(81)
sqrt(glengths)

round(3.14159)

?round
args(round)
example("round")

round(3.14159, digits = 2)

mean(glengths)

# Create the vector
test <- c(1, NA, 2, 3, NA, 4)

# Calculate the mean, ignoring NA values
mean_value <- mean(test, na.rm = TRUE)

# Print the result
print(mean_value)

# na.rm = TRUE → tells R to remove NA values before calculating the mean.

# Sort in descending order
sorted_glengths <- sort(glengths, decreasing = TRUE)

# Print the result
print(sorted_glengths)

square_it <- function(x) {
  square  <- x * x
  return(square)
}
square_it(5)


multiply_it <- function(x, y) {
  multiply <- x * y
  return(multiply)
}

multiply_it(4, 6)

?read.csv

metadata <- read.csv("mouse_exp_design.txt")

proj_summary  <- read.table("project-summary.txt", header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
head(proj_summary)


metadata

head(metadata)
str(metadata)
class(metadata)
summary(metadata)
tail(metadata)
length(metadata)
dim(metadata)
nrow(metadata)
ncol(metadata)
rownames(metadata)
colnames(metadata)


class(glengths)
class(metadata)

summary(proj_summary)

length(Samplegroup)

dim(proj_summary)

length(colnames(proj_summary))
length(proj_summary)
colnames(proj_summary)


# Excersice
temp_conv  <- function(x) {
  temp_c  <- (x-32) * 5/9
  temp_k  <- temp_c + 273.15
  return(temp_k)
}

temp_conv(70)

round(temp_conv(70))


age <- c(15, 22, 45, 52, 73, 81)
age[5]
age[-5]

age[c(3, 5,6)]

idx <- c(3, 5, 6)
age[idx]

age[1:4]

alphabets <- c("C", "D", "X", "L", "F")
edx <- c(1, 2, 5)
alphabets[edx]

alphabets[-3]
alphabets[5:1]

age[edx]

age
age > 50

age > 50 | age < 18
age

age[age > 50 | age < 18]

idx <- age > 50 | age < 18
age[idx]


which(age > 50 | age < 18)
age[which(age > 50 | age < 18)]

idx_num <- which(age > 50 | age < 18)
age[idx_num]

expression[expression == "high"]

idx <- expression =="high"
expression[idx]
idx

samplegroup[samplegroup == "KO"]

expression
str(expression)

samplegroup <- factor(samplegroup, levels = c("KO", "CTL", "OE"))

levels(samplegroup)


install.packages("ggplots2")
install.packages("ggplot2")

install.packages("BiocManager")

BiocManager::install("ggplot2")

install.packages("tidyverse")

# Extract balue 'Wt'
metadata[1, 1]

# Extract value '1'
metadata[1, 3]

# Extract third column
metadata[ , 3]

# Extract third column as a data frame
metadata[ , 3, drop = FALSE]

#Dataframe containing first two columns
metadata[ , 1:2]

# Data frame containing first, third and sixth rows
metadata[c(1,3,6), ]

# Extract the celltype column for the first three samples
metadata[c("sample1", "sample2", "sample3") , "celltype"]

#Extract the genotype column
metadata$genotype

metadata$genotype[1:5]

metadata[c("sample2", "sample8"), c("genotype", "replicate")]

metadata[c(4,9), "replicate"]

metadata$replicate

metadata$celltype == "typeA"

logical_idx <- metadata$celltype == "typeA"
logical_idx

metadata[logical_idx, ]

which(metadata$celltype == "typeA")

metadata[idx, ]
which(metadata$replicate > 1)

idx <- which(metadata$replicate > 1)
metadata[idx, ]

metadata[which(metadata$replicate >1), ]

sub_meta <- metadata[which(metadata$replicate >1), ]

Ko_meta <- metadata[which(metadata$genotype == "KO"), ]
Ko_meta

list1[[2]]

comp2 <- list1[[2]]
class(comp2)

list1[[1]]

list1[[1]][1]

random <- list(metadata, age, list1, samplegroup, number)
random

random[[4]]

names(list1)

# name components of the list
names(list1) <- c("species", "df", "number")
names(list1)

list1$df

rpkm_data <- read.csv("counts.rpkm.txt")
head(rpkm_data)

ncol(rpkm_data)
nrow(metadata)

A <- c(1, 3, 5, 7, 9, 11)
B <- c(2, 4, 6, 8, 1, 5)

# test to see if each of the elememts of A is in B
A %in% B

intersection <- A %in% B
intersection

A[intersection]

any(A %in% B)

all(A %in% B)

B[B %in% A]

A <- c(10, 20, 30, 40, 50)
B <- c(50, 40, 30, 20, 10)

# test to see if each element of A is in B
A %in% B

# test to see if each element of A is in the same position in B
A == B
# use all () to check if they are a perfect match
all(A == B)

x <- rownames(metadata)
y <- colnames(rpkm_data)

all(x %in% y)

all(rownames(metadata) %in% colnames(rpkm_data))

x == y
all(x == y)

important_genes <- c("ENSMUSG00000083700", "ENSMUSG00000080990", "ENSMUSG00000065619", "ENSMUSG00000047945", "ENSMUSG00000081010", "ENSMUSG00000030970")

rpkm_data[rpkm_data %in% important_genes]

rpkm_data[rownames(rpkm_data) %in% important_genes, ]

rpkm_data[important_genes, ]

teaching_team <- c("Jihe", "Mary", "Meeta", "Radhika", "Will", "Emma", "Heather", "Elizabeth", "Noor", "Upen")

# Extracting values  from a vector
teaching_team[c(2, 4)]

teaching_team

# Extracting all values and reordering them 
teaching_team[c(5, 4, 10, 6, 9, 2, 8, 1, 7, 3)]

# Saving the results to a variable
reorder_teach <- teaching_team[c(5, 4, 10, 6, 9, 2, 8, 1, 7, 3)]

first <- c("A", "B", "C", "D", "E")
second <- c("B", "D", "E", "A", "C")
second[c(4, 1, 5, 2, 3)]

match(first, second)

# Saving indices for how to reorder 'second' to match 'first'
reorder_idx <- match(first, second)

second[reorder_idx]

#Reordering and saving the output to a variable 
second_reordered <- second[reorder_idx]

first <- c("A", "B", "C", "D", "E")
second <- c("D", "B", "A")

match(first, second)

second[match(first, second)]

rownames(metadata)
colnames(rpkm_data)

genomic_idx <- match(rownames(metadata), colnames(rpkm_data))

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame 
rpkm_ordered <- rpkm_data[ , genomic_idx]

#View the reordered counts
View(rpkm_ordered)

all(rownames(metadata) == colnames(rpkm_ordered))

subset_rpkm <- rpkm_ordered[, !(colnames(rpkm_ordered) %in% c("sample2", "sample9"))]

View(subset_rpkm)

metadata_idx <- match(colnames(subset_rpkm), rownames(metadata))
subset_metadata <- metadata[metadata_idx, ]

mean(rpkm_ordered$sample1)

library(purrr) # Load the purrr
samplemeans <- map_dbl(rpkm_ordered, mean)

# Named vectors have a name assigned to each element instead of just referring to them as indices ([1], [2] and so on)
samplemeans

# Check length of the vector before adding it to the data frame
length(samplemeans)

# Create a numeric vector with ages. Note that there are 12 elements here 
age_in_days <- c(40, 32, 38, 35, 41, 32, 34, 26, 28, 28, 30, 32)

# Add the new vector as the last column to the new_metadata dataframe
new_metadata <- data.frame(metadata, samplemeans, age_in_days)

# Take a look at the new_metadata object
View(new_metadata)

animals <- read.csv("animals.txt")

is.data.frame(animals)

nrow(animals)
ncol(animals)

animals$speed == 40
animals[, "speed"]
animals[animals$speed == 40]
animals$speed[animals$speed == 40]

animals$color[animals$color == "Tan"]
animals$color == "Tan"
animals[animals$color == "Tan", ]

animals[which(animals$speed > 50), ]
subset(animals, speed > 50, select = color)

animals$color[animals$color == "Grey"] <- "Grey"

animal_list <- list(speed = animals$speed, color = animals$color)
animal_list
names(animal_list) <- c("speed", "color")
animal_list

ctrl_samples <- data.frame(row.names = c("sample3", "sample10", "sample8", "sample4", "sample15"), date = c("01/13/2018", "03/15/2018", "01/13/2018", "09/20/2018","03/15/2018"))

sum(rownames(ctrl_samples) %in% (proj_summary))

proj_summary_ctrl <- proj_summary[rownames(proj_summary) %in% rownames(ctrl_samples), ]
proj_summary_ctrl

matched_rows <- match(rownames(proj_summary_ctrl), rownames(ctrl_samples))
View(matched_rows)

ctrl_batch_info <- ctrl_samples[matched_rows, ]

proj_summary_ctrl <- cbind(proj_summary_ctrl, batch = ctrl_samples[matched_rows, "date"])
View(proj_summary_ctrl)

proj_summary_noctl <- proj_summary[proj_summary$treatment %in% c("high", "low"), ]

library(purrr)
proj_summary_noctl <- proj_summary_noctl[, map_lgl(proj_summary_noctl, is.numeric)]

## Load the new_metadata data frame into your environment from a .RData
load("new_metadata.RData")
# this data frame should have 12 rows and 5 columns
View(new_metadata)

library(ggplot2)

ggplot(new_metadata)

ggplot(new_metadata) + geom_point()

ggplot(new_metadata) + geom_point(aes(x = age_in_days, y = samplemeans))

ggplot(new_metadata) + geom_point(aes(x = age_in_days, y = samplemeans, colour = genotype))

ggplot(new_metadata) + geom_point(aes(x = age_in_days, y = samplemeans, colour = genotype, shape = celltype))

ggplot(new_metadata) + geom_point(aes(x = age_in_days, y = samplemeans, colour = genotype, shape = celltype), size = 2.25)

ggplot(new_metadata) + geom_point(aes(x = age_in_days, y = samplemeans, colour = genotype, shape = celltype), size = 3.0) + theme_bw()

ggplot(new_metadata) + geom_point(aes(x = age_in_days, y = samplemeans, color = genotype, shape = celltype), size = 2.25) + theme_bw() + theme(axis.title = element_text(size = rel(1.5)))


ggplot(new_metadata) + geom_point(aes(x = age_in_days, y = samplemeans, color = genotype, shape = celltype), size = 2.25) + 
theme_bw() + 
theme(axis.title = element_text(size = rel(1.5))) + 
xlab("Age (days)") +
ylab("Mean expression") +
ggtitle("Gene Expression vs Mouse Age") +
theme(plot.title = element_text(hjust = 0.5))

personal_theme <- function(){
  theme_bw() +
    theme(axis.title = element_text(size = rel(1.5))) +
    theme(plot.title = element_text(size = rel(1.5), hjust = 0.5))
}


ggplot(new_metadata) + 
  geom_point(aes(x = age_in_days, y = samplemeans, color = genotype, shape = celltype), size = rel(3.0)) +
  xlab("Age (days)") +
  ylab("Mean expression") +
  ggtitle("Expression with Age") +
  personal_theme()

ggplot(new_metadata) +
  geom_boxplot(aes(x =genotype, y = samplemeans, fill = celltype, )) +
  ggtitle("Genotype differences in average gene expression") +
  xlab("Genotype") +
  ylab("Mean expression") +
  theme_bw() +
  theme(
    axis.title = element_text(size = rel(1.25)),
    plot.title = element_text(size = rel(1.5), hjust = 0.5)
  )
new_metadata$genotype <- factor(new_metadata$genotype, levels = c("Wt", "KO"))

library(ggplot2)

ggplot(new_metadata) +
  geom_boxplot(aes(x = genotype, y = samplemeans, fill = celltype)) +
  scale_fill_manual(values = c("purple", "orange")) +
ggtitle("Genotype difference in average gene expression") +
  xlab("Genotype") +
  ylab("Mean expression") +
  theme_bw() +
  theme(
    axis.title = element_text(size = rel(1.25)),
    plot.title = element_text(size = rel(1.5), hjust = 0.5)
  )

ggplot(new_metadata) +
    geom_boxplot(aes(x = genotype, y = samplemeans, fill = celltype)) +
    scale_fill_manual(values = c("#800080", "#FFA500")) +  # Hex codes for purple and orange
    ggtitle("Genotype difference in average gene expression") +
    xlab("Genotype") +
    ylab("Mean expression") +
    theme_bw() +
    theme(
      axis.title = element_text(size = rel(1.25)),
      plot.title = element_text(size = rel(1.5), hjust = 0.5)
    )

# Save a data frame to file
write.csv(sub_meta, file = "subset_meta.txt")

?write.csv

# Save a vector to file
write(glengths, file = "genome_lengths.txt")

# Save a vector to file as a single column
write(glengths, file = "genome_lengths.txt", ncolumns = 1)

ggplot(new_metadata) + 
  geom_point(aes(x = age_in_days, y = samplemeans, color = genotype, shape = celltype), size = rel(3.0))
ggsave("figures/scatterplot.pdf")


# Create vector of work days
work_days <- c(Monday, Tuesday, Wednesday, Thursday, Friday)

# Create a function to round the output of the sum function 
round_the_sum <- function(x){
  return(round(sum(x)))
}

# Create a function to add together three numbers
add_numbers <- function(x,y,z){
  sum(x,y,z)
}
add_numbers(5,9)

Error: package or namespace load failed for 'Seurat' in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]): there is no package called 'multtest'


library(tidyverse)

## A single command
sqrt(83)

## Base R method of running more than one command
round(sqrt(83), digits = 2)

## Running more than one command with piping 
sqrt(83) %>% round(digits = 2)

random_numbers <- c(81, 90, 65, 43, 71, 29)

mean(random_numbers)
round(mean(random_numbers), digits = 2)

##Running more than one command with piping
mean(random_numbers) %>% round(digits = 2)

library(tidyverse)

# 1. Read the g:Profiler results file
gprofiler_data <- read_tsv("gprofiler_results_Mov10oe.txt")

# 2. Filter for GO Biological Process (BP) terms
bp_data <- gprofiler_data %>%
  filter(domain == "BP")

# 3. Prepare and rename relevant columns
go_bp_clean <- bp_data %>%
  select(
    Term = term.name,
    GO_ID = term.id,
    PValue = p.value,
    GeneCount = overlap.size,
    TermSize = term.size,
    QuerySize = query.size
  ) %>%
  mutate(
    GeneRatio = GeneCount / TermSize,
    logP = -log10(PValue)
  )

# 4. Take the top 30 most significant GO terms
top_go <- go_bp_clean %>%
  arrange(PValue) %>%
  slice(1:30)

# 5. Plot: Dotplot of top GO BP terms
ggplot(top_go, aes(x = GeneRatio, y = reorder(Term, GeneRatio))) +
  geom_point(aes(size = GeneCount, color = logP)) +
  scale_color_gradient(low = "yellow", high = "red", name = "-log10(p-value)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "Top 30 Enriched GO:BP Terms",
    x = "Gene Ratio",
    y = "GO Term"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10)
  )

# Read in the fucntional analysis results
functional_GO_results <- read_delim(file = "gprofiler_results_Mov10oe.txt", delim = "\t")

functional_GO_results


# Return only GO biological processes
bp_oe <- functional_GO_results %>%
  filter(domain =="BP")

view(bp_oe)

library(dplyr)

bp_oe <- bp_oe %>%
  filter(relative.depth > 4)

# Selecting columns to keep
bp_oe <- bp_oe  %>%
  select(term.id, term.name, p.value, query.size, term.size, overlap.size, intersection)

view(bp_oe)

# Provide better names for columns
bp_oe <- bp_oe %>%
  dplyr::rename(GO_id = term.id,
                GO_term = term.name)

# Create gene ratio column based on other columns in dataset
bp_oe <- bp_oe  %>%
  mutate(gene_ratio = overlap.size / query.size)

bp_oe <- bp_oe %>%
  mutate(term_percent = (overlap.size / term.size) * 100)


# Load necessary package
library(tibble)

# Convert data frame to tibble and add row names as a column
animals_tb <- rownames_to_column(animals, var = "animal_names")

# Load ggplot2
library(ggplot2)

# Create scatterplot
ggplot(animals_tb, aes(x = animal_names, y = speed)) +
  geom_point(color = "blue", size = 3) +               # Scatterplot with blue points
  theme_minimal() +                                    # Clean minimal theme
  labs(
    title = "Animal Speeds",
    x = "Animal",
    y = "Speed (km/h)"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1), # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5)             # Center title
  )

  )

