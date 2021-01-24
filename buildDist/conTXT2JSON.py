#!/usr/bin/python
# import packages to read a csv file and write a new json file

import csv, json

csvFilePath = "Gen118mod.csv"
jsonFilePath = "Gen118.json"

# Read CSV file and add to data
data = {}
with open(csvFilePath) as csvFile:
    csvReader = csv.DictReader(csvFile)
    for rows in csvReader:
        print(rows)
        ConnNode = rows["ConnNode"]
        data[ConnNode] = rows

# Create new JSON file and add data to it
with open(jsonFilePath, "w") as jsonFile:
    # make it more readable and pretty
    jsonFile.write(json.dumps(data, indent=4))

