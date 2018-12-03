Implementation of alignment-free sequence comparison methods - motivation and results
=========================================================================
<b>Motivation:</b>b> Alignment-based methods for sequence comparison proved to be fundamental to clarifying the function of thousands of various protein or gene families. However, processes such as recombination, exon/domain shuffling or the existence of low sequence similarity families make these methods ineffective. As a result, a new sequence comparison methods were developed, which provide alternatives over alignment-based approaches. Implementations of alignment-free methods are still rare, but if they exist: (i) they are not included in a single program or website, (ii) they use external programs, and (iii) their service is extremely difficult and non-intuitive for the average user.
</br></br>
<b>Results:</b>b> Fourteen alignment-free methods were developed as a result of this paper. All of them are part of "Alfree" meta-server, which is being constantly developed in Laboratory of Computational Biology, Adam Mickiewicz University, Poznań, Poland. “Alfree” facilitates the calculation and visualization of the evolutionary relationships between sequences. Developed software for alignment-free sequence comparisons was also made available as scripts and library written in the Python programming language, at: https://github.com/krzyszto9/alignment-free. Furthermore, within this paper there was performed the evaluation of predictions of developed software’s quality. It was performed by comparing the predictions to a reference set of protein families which is deposited in the SCOPe database, and by analyzing the ROC curves. Euclidean distance and the square Euclidean distance had the highest sensitivity and specificity of predictions (AUC value= 0.666). Whereas the method that had the lowest values of both parameters was normalized compression distance (AUC = 0.555). Compared to the alignment-free methods, the Smith-Waterman algorithm had the best AUC values for three out of four analyzed groups of sequences. Running time analysis showed that methods implemented in this paper were over thousand times faster than the Smith-Waterman algorithm.
</br>

Implemented distances/methods
===============================
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

Full text
============
Full text of my Bachelor thesis paper "Implementation of alignment-free sequence comparison methods" is available (in Polish) under this address:</br>
<a href="https://github.com/krzyszto9/alignment-free/blob/master/BSc_aligment_free.pdf">Krzysztof Udycz, Implementation of alignment-free sequence comparison methods</a>
</br>

Help
====
The --help/-h shows a description of what each selected script does.
</br>
