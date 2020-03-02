# Concatenate openQCD files
---

Script to concatenate openQCD files into a single file. The script generates
a file for each particle with a given gamma structure and external momenta.

In order to run the file, you shall use the following command,

```bash
python ./concatenate_files.py input_dir output_dir ( N_gammas )
```
where the `input_dir` is the one produced by openqcd_hadspec, normally 
named `dat`. Besides, `output_dir` is optional and arbitrary and sets the
directory where the data will be stored. The last parameter sets the number 
of gamma structures, if not set, it uses `16` gamma structues. It is useful 
for non-diagonal correlation functions.

After analysing the data, it is scattered inside `output_dir`, classified by
particle `[meson, baryon]`, then gamma structure or baryon structure. Inside
each particle, it is classified by value of momenta `p2Max`.


