function write_tsv(tab, tsv_file)
%WRITE_TSV   Write a table into a BIDS tab-separated value file.
%
%  write_tsv(tab, tsv_file)

writetable(tab, tsv_file, 'FileType', 'text', 'delimiter', '\t');
