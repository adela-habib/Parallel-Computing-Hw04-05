import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np



arr = np.zeros((64,64))


line_index = 0

for line in open("0.25HMtest.csv"):

	if line_index > 63:
		break

	nums = line.split(",")

	# get rid of null bytes
	if(nums[0].startswith("\x00")):
		nums[0] = nums[0].replace("\x00", "")
		

	num_index = 0
	for val in nums:

		if val == "\n":
			break
		try:
			int_val = int(val)
		except:
			print("ERROR: ", val)
			print("line: ", line_index)

			print("num: ", num_index)


		arr[line_index][num_index] = int_val/float(256)	

		num_index +=1




	line_index += 1 

print("Done parsing, creating heatmap")

ax = sns.heatmap(arr, xticklabels=False, yticklabels=False)
plt.show()





