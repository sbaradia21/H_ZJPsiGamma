#!/bin/bash
new_words=("GluGlu" "VBF" "WH" "ZH" "BBH" "TTH") # array of new words
for new_word in "${new_words[@]}"; do
  new_filename="make_config_2018_${new_word}.sh" # create a new filename for each copy
  cp make_config_2018.sh "${new_filename}" # copy the original file
  sed -i "s/orange/${new_word}/g" "${new_filename}" # replace the word "apple" with the new word
done
