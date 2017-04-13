if [ $# -eq 0 ]
  then
    echo "Please supply the name of the database as the first argument. Either small-data or full-data"
else
  export SNORKELDB="postgres:///snorkel-$1"
  export AGP_DATA_SIZE="$1"
fi
