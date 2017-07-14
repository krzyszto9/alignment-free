Implementation of alignment-free sequence comparison methods - motivation
=========================================================================
Alignment-based methods for sequence comparison proved to be fundamental to clarifying the function of thousands of various protein or gene families. However, processes such as recombination, exon/domain shuffling or the existence of low sequence similarity families make these methods ineffective. As a result, a new sequence comparison methods were developed, which provide alternatives over alignment-based approaches. Implementations of alignment-free methods are still rare, but if they exist: (i) they are not included in a single program or website, (ii) they use external programs, and (iii) their service is extremely difficult and non-intuitive for the average user.
</br>

Implemented distances/methods
===================
Measures based on word frequencies (file ltuple.py):
  - eucliudean distance;
  - squared eucliudean distance;
  - standardized eucliudean distance;
  - standardized eucliudean distance with different resolutions;
  - weighted eucliudean distance;
  - weighted eucliudean distance with different resolutions;
  - minkowski distance;
  - Pearson productmoment correlation coefficient;
  - Kullback–Leibler divergence;
  - cosine distance;
  - evolutionary distance.
  
Information theory-based measures (file infotheory.py):
  - Universal Sequence Maps;
  - normalized compression distance.

More informations about alignment-free methods:</br>
<a href="https://github.com/krzyszto9/alignment-free/blob/master/Vinga%2C%20Almeida%20-%202003%20-%20Alignment-free%20sequence%20comparison--a%20review.pdf">Vinga, S.; Almeida, J. (Mar 1, 2003) Alignment-free sequence comparison-a review.
Bioinformatics 19 (4): 513–23</a>
</br>

Help
====
The --help/-h shows a description of what each selected script does.
