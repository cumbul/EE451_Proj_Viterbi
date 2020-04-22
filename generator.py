# data generator
import numpy as np


# num1 -- number of states
# num2 -- number of observations inside state
def generate(num1, num2):
	# generate states
	trans_prob = dict(dict())
	for i in range(num1):
		name_from = "S" + str(i)
		probs = np.random.random(num1)
		probs /= probs.sum()
		for j in range(num1):
			name_to = "S" + str(j)
			trans_prob[name_from] = {name_to:probs[j]}

	# generate observations
	observation_prob = dict(dict())
	for i in range(num1):
		name_from = "S" + str(i)
		probs = np.random.random(num2)
		probs /= probs.sum()
		for j in range(num2):
			name_to = "O" + str(j)
			observation_prob[name_from] = {name_to:probs[j]}

	# generate initial probability
	init_prob = dict()
	probs = np.random.random(num1)
	probs /= probs.sum()
	for i in range(num1):
		name = "S" + str(i)
		init_prob[name] = probs[i]

	file1 = open("data_transition.txt","w+") 
	file1.write(str(trans_prob))
	file1.close()

	file2 = open("data_observation.txt","w+") 
	file2.write(str(observation_prob))
	file2.close()

	file3 = open("data_init.txt","w+") 
	file3.write(str(init_prob))
	file3.close()

generate(100,50)