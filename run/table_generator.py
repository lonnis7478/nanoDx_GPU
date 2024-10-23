import pandas as pd
import os

def sort(values):
    return values.str.extract('(\d+)$')[0].astype(int)


#path="/home/sander/Documents/sturgeon/result/"

#table = pd.DataFrame(columns = ["Sample", "Modkit", "Bed", "Prediction", "Total"])

#for folder in os.listdir(path):

#	file = open(path+folder+"/benchmark.txt")
#	data = file.read().split("\\n")[1].split("\\t")
#	new_row = {"Sample": folder, "Modkit":int(data[0]), "Bed":int(data[1]), "Prediction":int(data[2]), "Total":int(data[3][:-1])}
#	table = pd.concat([table, pd.DataFrame([new_row])], ignore_index=True)
#	file.close()
	
#table[["Modkit", "Bed", "Prediction", "Total"]] /= 1000
#table = table.sort_values(by=["Sample"],key=sort)
#print(table)
#table.to_csv("sturgeon_time.csv")


path="../results/classification/"

table = pd.DataFrame(columns = ["Sample", "cpgs", "Score", "Prediction"])

for file in os.listdir(path):
	if ".txt" in file:
		sample = file.split("-")[0]

		votes = open(path+file)
		data = votes.read().split("\n")
		cpgs = int(data[0].split(":")[1])
		pred = data[1].split(":")[1]
		score = data[2].split(":")[1]

		new_row = {"Sample": sample, "cpgs": cpgs, "Score": score, "Prediction": pred}
		table = pd.concat([table, pd.DataFrame([new_row])], ignore_index=True)


		votes.close()
		
table = table.sort_values(by=["Sample"],key=sort)
print(table)
print("(",table[table["cpgs"] < 3000]["Sample"].values,")")
