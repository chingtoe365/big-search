# big-search
Case Study in Python searching datasets larger than RAM.

{{ Describe your implementation here}}


- First seqarate sequence from different chromosome into different txt file
- For each chromosome, split the sequence into chunks, create suffix array for each trunk (partition) using an efficient algorithm designed by Mark Mazumder with complexity of O(n)
	### Several methods have been tried to applied here
	### the SA-IS algorithm is able to create a suffix array for the whole chromosome without causing memory issues, but limitted by my own machine the huge suffix array file cannot be read afterwards in the profile search so I skip this method
	### a self-innovated algorithm was designed here as well, with less memory consumption, but slower in speed
- For each character in the pattern, for each trunk in each chromosome, perform a binary search to find the intervals which contain the next matching character, and then alternatively increment the mismatch counter for the candidate match
	### recursively using the same suffix array calculated beforehand for this trunk (partion), with just some minor changes to fit with the order of characters
- Remove candidate matches that exceed k mismatches.
- Stop when the pattern is fully traversed

! - There exist differences between bowtie read sequences and the sequence I have got with the exact location (bowties' might have been modified).
