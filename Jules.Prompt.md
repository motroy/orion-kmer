- Read the SPECS.md, generate a tasks list in TODO.md for steps required to produce the code for the specified tool. Then follow the tasks in TODO.md and generate the code for the tool, checking off the completed tasks as you go
- add function classify - user gives genome or reads and orion-kmer compares to database (or multiple databases) generating output of proportion of matches in the database and their coverage (depth and breadth)
  - that's ok, should also include a breadth of coverage that looks at the Number of unique k-mers from the input that were found in that specific match in the database (instead of all the database)
  - per individual genome/contig within a database if a database can represent multiple references
- Add a progress bar when running a function, and print execution time to stdout when run completes. all functions in orion-kmer. Also print to stdout the max ram usage
- For orion-kmer classify function, add option minimum coverage to show in output (i.e. don't show db reference genomes that had no kmers matching or zero coverage). Also add an optional output tsv that shows the results summary.
