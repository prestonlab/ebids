function tab = read_tsv(tsv_file)
%READ_TSV   Read a BIDS tab-separated value file into a table.
%
%  tab = read_tsv(tsv_file)

tab = readtable(tsv_file, 'FileType', 'text', 'delimiter', '\t', ...
                'TreatAsEmpty', 'n/a');
