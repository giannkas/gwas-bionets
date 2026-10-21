# samfam, K, I

# Generate a pseudorandom byte stream. Using openssl command and AES-256-CTR cipher.
get_seeded_random()
{
    seed="\$1"
    openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt \
    </dev/zero 2>/dev/null
}

# Compute number of lines (non-empty lines) in the samfam file.
total_samples=\$(grep -c '^' $samfam)

# Number of samples to be included in each subsample.
n_samples=\$(expr \$total_samples - \$total_samples / $K)

# We extract different basename (without extension) depending on the PLINK version.
if [ $vplink -eq 1 ]; then
    nosuffix=\$(basename -s .fam $samfam)
else
    nosuffix=\$(basename -s .psam $samfam)
fi

# We add a suffix split if K is greater than 1, otherwise we are not splitting.
if [ $K -gt 1 ]; then
    split_suffix="_split_${I}"
    n_samples=\$(expr \$total_samples - \$total_samples / $K)
else
    split_suffix=""
    n_samples=\$(expr \$total_samples)
fi

# We create the subsample files by shuffling lines and according to n_samples number.
# It creates different file formats whether PLINK is 1 or 2.
if [ $vplink -eq 1 ]; then
    cut -f1,2 -d' ' $samfam | shuf --random-source=<(get_seeded_random ${I}) | tail -n \$n_samples > \${nosuffix}\${split_suffix}.fam
else
    cut -f1 $samfam | shuf --random-source=<(get_seeded_random ${I}) | tail -n \$n_samples > \${nosuffix}\${split_suffix}.psam
fi