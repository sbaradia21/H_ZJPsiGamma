#!/bin/bash
new_words=("sigmaup" "sigmadown" "scaleup" "scaledown" "muup" "mudown" "nominal") # array of new words
for new_word in "${new_words[@]}"; do
  new_filename="HtoMuMuG_ZH_2018_${new_word}.job" # create a new filename for each copy
  cp HtoMuMuG_ZH_2018.job "${new_filename}" # copy the original file
  sed -i "s/apple/${new_word}/g" "${new_filename}" # replace the word "apple" with the new word
done
