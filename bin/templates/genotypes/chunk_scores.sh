# SCORES, K, I, chunk_suffix

# generate pseudorandom byte stream
get_seeded_random()
{
    seed="\$1"
    openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt \
    </dev/zero 2>/dev/null
}

# count the number of non-empty lines minus 1 (the header)
total_samples=\$(( \$(grep -c '^' ${SCORES}) - 1 ))

# removing file extension
nosuffix=\$(basename -s .tsv ${SCORES})

# changing chunk suffix depending on K value
# computing the number of samples to save in each subsample file  (n_samples)
if [ $K -gt 1 ]; then
    chunk_suffix="_chunk_${I}"
    n_samples=\$(expr \$total_samples - \$total_samples / $K)
else
    chunk_suffix=""
    n_samples=\$(expr \$total_samples)
fi

# extracting the header given that we need to include in each subsample file
header=\$(head -n 1 ${SCORES})

# extract the random subsamples without the header 
shuf --random-source=<(get_seeded_random ${I}) ${SCORES} | tail -n +2 | head -n \$n_samples >\${nosuffix}\${chunk_suffix}.tsv

# prepend the header to the subsamples file
sed -i "1i\${header}" \${nosuffix}\${chunk_suffix}.tsv