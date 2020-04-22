# data generator
import numpy as np
import json

# num1 -- number of states
# num2 -- number of observations inside state
def generate(num1, num2):
	# generate states
	trans_prob = dict(dict())
	for i in range(num1):
		name_from = "State" + str(i)
		probs = np.random.random(num1)
		probs /= probs.sum()
		trans_prob[name_from] = {}
		for j in range(num1):
			name_to = "State" + str(j)
			trans_prob[name_from].update({name_to:probs[j]}) 

	# generate observations
	observation_prob = dict(dict())
	for i in range(num1):
		name_from = "State" + str(i)
		probs = np.random.random(num2)
		probs /= probs.sum()
		observation_prob[name_from] = {}
		for j in range(num2):
			name_to = "Obs" + str(j)
			observation_prob[name_from].update({name_to:probs[j]})

	# generate initial probability
	init_prob = dict()
	probs = np.random.random(num1)
	probs /= probs.sum()
	for i in range(num1):
		name = "State" + str(i)
		init_prob[name] = probs[i]

	str_trans  = json.dumps(trans_prob, indent=4, sort_keys=False)
	str_observation  = json.dumps(observation_prob, indent=4, sort_keys=False)
	str_init  = json.dumps(init_prob, indent=4, sort_keys=False)

	file1 = open("data_transition.txt","w+") 
	file1.write(str_trans)
	file1.close()

	file2 = open("data_observation.txt","w+") 
	file2.write(str_observation)
	file2.close()

	file3 = open("data_init.txt","w+") 
	file3.write(str_init)
	file3.close()

generate(100,50)