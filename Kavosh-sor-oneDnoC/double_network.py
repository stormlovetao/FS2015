##This program change a network from direct graph to undirect one

filein = open(r'./network/V4_net_exp.txt','r')
fileout = open(r'./network/V4_net_exp_double.txt','w')
firstline = filein.readline()
if len(firstline.strip()) == 0:
	print "there is no nodes there! Exit!"
	exit()
fileout.write(firstline)
while True:
	line = filein.readline()
	if len(line) == 0:
		break
	#print line.strip()
	linesplit = line.strip().split('\t')
	reverse_line = linesplit[1] + '\t' + linesplit[0] + '\n'
	fileout.write(line)
	fileout.write(reverse_line)
