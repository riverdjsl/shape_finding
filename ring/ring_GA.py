import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def trans_comma_str(string):
	'''to transform comma strings into lists of numbers'''
	item = []
	lst = []
	for i in string[:-1]:
		if i != ",":
			item.append(i)
		else:
			lst.append(''.join(item))
			item = []
			continue
	return lst


class Steel345():
	def __init__(s, t):
		s.t = t

	def fy(s):
		if s.t <= 16:
			return 305
		elif 16 < s.t <= 40:
			return 295
		elif 40 < s.t <= 63:
			return 290
		elif 63 < s.t <= 80:
			return 280
		elif 80 < s.t <= 100:
			return 270
		else:
			return None

	def E(s):
		return 206e3


def mtgenerater(ts):
	return [Steel345(ts[0]), Steel345(ts[1])]


class Pipe():
	def __init__(s, D, t):
		s.D = D
		s.t = t
		s.r = D/2

	def A(s):
		return np.pi*s.D*s.t

	def I(s):
		return np.pi*(s.D/2)**3*s.t

	def W(s):
		return np.pi*(s.D/2)**2*s.t

	def i(s):
		return np.sqrt((np.pi*(s.D/2)**3*s.t)/(np.pi*s.D*s.t))


def pipegenerater(ts, factor):
	return [Pipe(ts[0]*factor, ts[0]), Pipe(ts[1]*factor, ts[1])]


def fi(lambdai, f, E):
	alpha1, alpha2, alpha3 = 0.41, 0.986, 0.152
	lambdan = lambdai*np.sqrt(f/E)/np.pi
	if 0 < lambdan <= 0.215:
		return 1-alpha1*lambdan**2
	elif lambdan > 0.215:
		return ((alpha2+alpha3*lambdan+lambdan**2)-\
		np.sqrt((alpha2+alpha3*lambdan+lambdan**2)**2-4*lambdan**2))\
		/(2*lambdan**2)
	

def check2(ratio, N, l, u3, u2, E, f, W, I, A, i, r, t):

	lambdax = u3*l/i
	lambday = u2*l/i
	
	if lambdax < 0.1:
		Nex = 99999
		fix = 1.0
		fiy = 1.0
	else:
		Nex = np.pi**2*E*I/(1.1*lambdax**2)
		fix = fi(lambdax, f, E)
		fiy = fi(lambday, f, E)
	
	part1 = N/(min(fix, fiy)*A*f)
	part4 = N/(A*f)
	condition1 = abs(part4)
	condition2 = abs(part1)
	
	if condition1 <= ratio and condition2 <= 1 and max(lambdax, lambday) <= 120*np.sqrt(235/345):
		return True
	else:
		return False


def strainenergy(N, l, u3, u2, E, f, W, I, A, i, r, t):
	'''pipe only'''

	part1 = (N)**2*l/(2*E*A)
	Total_energy = part1
	Volume = A*l

	return Total_energy, Volume


def vars(pp, mt, mu):
	return {'u3': mu[0], 'u2': mu[1], 'E': mt.E(), 'f': mt.fy(), \
			'W': pp.W(), 'I': pp.I(), 'A': pp.A(), 'i': pp.i(), \
			'r': pp.r, 't': pp.t}


def pointlist2(x1, x2, y2):
	return [(x1, y2), (x2, y2), (x2, 0)]


def chromosome_generator(rdn):
	return np.binary_repr(int(rdn*511), width=9)


def DEC(bi):

	result = []
	for i, j in enumerate(bi):
		result.append(int(j)*2**i)

	return sum(result)


def individual_generator(rdns):

	homo = []
	for i in rdns:
		homo.append(chromosome_generator(i))

	return tuple(homo)


def chromosome_translator(individual, tmax, start_elev, A, B, PC):

	def coord(chromosome, a, b):
		return a+(b-a)*DEC(chromosome)/511

	yi = coord(individual[0], start_elev, A[1])
	i = (PC[0], yi)

	xj = coord(individual[1], A[0], PC[0])
	yj = coord(individual[2], i[1], A[1])
	j = (xj, yj)

	t1 = coord(individual[3], 5, tmax[0])
	t2 = 2*t1

	return ([A, j], [B, j], [j, i], [i, PC], [A, B]), (t1, t2), (i, j)


#sfitness points: the minimun value gets the highest point.	
def point1(v, rou, TF, shears):
	'''Fitness based on total volumes!'''
	if TF is True:
		return 5-v/1e8
	else:
		return 0


def point2(v, rou, TF, shears):
	'''Fitness based on beam shear similarity.'''
	if TF is True:
		if abs(shears[0][0]-shears[1][0])<1:
			return 5-v/1e8
		else:
			return 0
	else:
		return 0


def funcchoose(funclist):
	for i, j in enumerate(funclist):
		print('Method ', i, ':', j.__doc__)
	while True:
		try:
			fn = input("Now choose one:")
			for i, j in enumerate(funclist):
				if i == int(fn):
					return j, j.__doc__
		except Exception:
			continue

f = [point1, point2]
funcchosen = funcchoose(f)
point = funcchosen[0]


#selecting high points chromosomes
def select(population, points):

	mark_a = np.random.rand()*sum(points)
	mark_b = ((mark_a+0.5)%1)*sum(points)

	parents = []
	for mark in [mark_a, mark_b]:
		accu = 0
		for i, j in zip(population, points):
			accu += j
			if accu >= mark:
				parents.append(i)
				break

	return parents


#genatic exchange
def genatic_exchange(chromosome_f, chromosome_m):
	chromosome_length = min(len(chromosome_f), len(chromosome_m))
	exchange_point1 = int(chromosome_length*np.random.rand())
	exchange_point2 = int(chromosome_length*np.random.rand())
	f = list(chromosome_f)
	m = list(chromosome_m)
	m[exchange_point1: exchange_point2], f[exchange_point1: exchange_point2] = \
		f[exchange_point1: exchange_point2], m[exchange_point1: exchange_point2]
	return ''.join(m)


#mutation
def mutation(chromosome):
	chromosome_length = len(chromosome)
	mutation_point = int(chromosome_length*np.random.rand())
	chance = np.random.rand()
	if chance >= 0.3:
		ch = list(chromosome)
		for i, gene in enumerate(ch):
			if  i == mutation_point:
				ch[i] = str((1-int(gene)))
		return ''.join(ch)
	else:
		return chromosome


def population_mutation(population, possibility=0.3):
	if np.random.rand() >= 1-possibility:
		mu_person_index = int(np.random.rand()*len(population))
		mu_person = population.pop(mu_person_index)

		new_indi = []
		for chromosome in mu_person:
			new_indi.append(mutation(chromosome))
		population.append(tuple(new_indi))
	return population


#reproduction
def reproduction(parents, children_num=4):

	children = []
	for i in range(children_num):

		embryo = []
		for i, j in zip(parents[0], parents[1]):
			embryo.append(genatic_exchange(i, j))

		child = []
		for i in embryo:
			child.append(mutation(i))

		children.append(tuple(child))

	return children


def rawdata(filein):

	datain = []
	with open(filein) as f:
		for i in f.readlines():
			datain.append(trans_comma_str(i))

	dataindic = {}
	for j in datain:
		dicj = {}
		dicj[j[0]] = [float(k) for k in j[1:]]
		dataindic.update(dicj)


	mus = [dataindic['mub'], dataindic['mut']]
	plist = pointlist2(*dataindic['geometry'])
	forces = dataindic['forces']
	ratio = dataindic['ratio']
	tmax = dataindic['tmax']
	
	start_elev = dataindic['starting_elevation'][0]

	p_num = int(dataindic['population'][0])
	g_num = int(dataindic['generation'][0])

	def energy(individual, F):

		individual_values = chromosome_translator(individual, tmax, start_elev, *plist)
		candidate = individual_values[0]
		pipets = individual_values[1]
		positions = individual_values[2]

		material = mtgenerater(pipets)
		pipes = pipegenerater(pipets, dataindic['D_over_t'][0])

		vars1, vars2 = (vars(i, j, k) \
		for i, j, k in zip(pipes, material, mus))

		energyi = []
		volumei = []
		checks = []
		
		def geo(line):
			'''
			line is [(x1, y1), (x2, y2)]
			'''
			dx = abs(line[0][0]-line[1][0])
			dy = abs(line[0][1]-line[1][1])
			try:
				ang = np.arctan(dx/dy)
			except ZeroDivisionError:
				ang = np.pi/4
			l = np.sqrt(dx**2+dy**2)
			return ang, l
		
		angs = []
		ls = []	
		for i in candidate:
			angs.append(geo(i)[0])
			ls.append(geo(i)[1])
		
		def solver(angs, ls, q):
			'''
			angs is [alpha, beta, gamma]
			ls is [ljA, ljB, lij, loi, lAB]
			q is [q]
			return (NjA, NjB, Nij, Noi)
			'''
			Rcy = q[0]*ls[4]
			Rcx = Rcy*np.tan(angs[2])
			A = np.array([[np.cos(np.pi/2-angs[0]), -np.cos(np.pi/2-angs[1])], \
						[np.sin(np.pi/2-angs[0]), np.sin(np.pi/2-angs[1])]])
			B = np.array([Rcx, Rcy])
			C = np.linalg.solve(A, B)
			dic1 = {'jA': (C[0], ls[0]), 'jB': (C[1], ls[1]), 'ij': (Rcy/np.cos(angs[2]), ls[2])}
			dic2 = {'oi': (2*Rcy, ls[3])}
			shear1 = C[0]*np.cos(angs[0])
			shear2 = C[1]*np.cos(angs[1])
			Mb0 = q[0]*ls[4]**2/12
			deltaM = shear1*ls[4]-q[0]*ls[4]**2/2
			MA = Mb0+deltaM/2
			MB = Mb0-deltaM/2
			return dic1, dic2, [(shear1/1000, MA/1e6), (shear2/1000, MB/1e6)]
			
		dic1 = solver(angs, ls, forces)[0]
		dic2 = solver(angs, ls, forces)[1]
		shears = solver(angs, ls, forces)[2]
		
		for i, j in dic1.items():
			energyi.append(strainenergy(*j, **vars1)[0])
			volumei.append(strainenergy(*j, **vars1)[1])
			checks.append(check2(ratio[0], *j, **vars1))
		for i, j in dic2.items():
			energyi.append(strainenergy(*j, **vars2)[0])
			volumei.append(strainenergy(*j, **vars2)[1])
			checks.append(check2(ratio[0], *j, **vars2))

		energy = sum(energyi)
		volume = sum(volumei)
		density = energy/volume
		
		return volume, density, all(checks), shears

	population = []
	for i in range(p_num):
		population.append(individual_generator(np.random.rand(4)))

	generation = 0
#	energy_g = list(map(energy, population, [forces]*len(population)))

	volumemins = []
	while generation <= g_num:

		print('\n\n\nGeneration {}.\n'.format(generation))

		volume_g = list(map(energy, population, [forces]*len(population)))
		pts = [point(i, j, k, ii) for i, j, k, ii in volume_g]

		if generation > 0.6*g_num and abs(np.average(pts)-max(pts)) <= 0.005:
			print("Ending the run 'cause convergence almost achieved!")
			break

		print('\n\n'.join(list('  '.join(list('{:>5.2f}'.format(i) for i in pts)[j: j+5]) for j in range(0, len(pts), 5))))

		pama = select(population, pts)
		kids = reproduction(pama)
		v_kids = list(map(energy, kids, [forces]*len(kids)))
		pts_k = [point(i, j, k, ii) for i, j, k, ii in v_kids]

		for i in sorted(pts[:])[:4]:
			for j in pts_k[:]:
				if j >= i:
					population[pts.index(i)] = kids[pts_k.index(j)]
					pts_k.pop(pts_k.index(j))

		if abs(np.average(pts)-max(pts)) > 1:
			population = population_mutation(population)
		else:
			population = population_mutation(population, possibility = 0.05)

		volumemins.append(min(volume_g)[0])
		print('\nCurrent Volumn: {:.1f}.'.format(min(volume_g)[0]))
		individual_em = population[volume_g.index(min(volume_g))]
		shearforces = min(volume_g)[3]

		generation += 1

	return chromosome_translator(individual_em, tmax, start_elev, *plist), dataindic['D_over_t'][0], min(volume_g)[0], volumemins, generation, shearforces


rd = rawdata('ring_GA.txt')

data = []
for i in rd[0][0][:-1]:
	xs = []
	ys = []
	for j in i:
		xs.append(j[0])
		ys.append(j[1])
	data.append([xs, ys])
	
for i in rd[0][0][:-2]:
	xs = []
	ys = []
	for j in i:
		xs.append(-j[0])
		ys.append(j[1])
	data.append([xs, ys])

pipe_description = []
for i in rd[0][1]:
	pipe_description.append('D = {:.1f}, t = {:.1f}'.format(i*rd[1], i))

shear_forces_description = []
for i, j in zip(rd[5], ['A', 'B']):
	shear_forces_description.append('V{} = {:.1f}, M{} = {:.1f}'.format(j, i[0], j, i[1]))

with plt.xkcd():
	fig1 = plt.figure('fig1')
	gs = GridSpec(3, 3)

	ax1 = plt.subplot(plt.subplot(gs[0:-1, 0:-1]))
	ax1.axis('equal')

	ax4 = plt.subplot(plt.subplot(gs[0:-1, 2]))
	ax4.text(0.0, 0.6, 'The total volume:', fontsize=10)
	ax4.text(0.1, 0.5, int(rd[2]), fontsize=20, color='red')

	ax4.text(0.6, 0.6, 'The uper branches:', fontsize=10)
	ax4.text(0.6, 0.55, pipe_description[0], fontsize=10, color='blue')
	ax4.text(0.6, 0.4, 'The trunk:', fontsize=10)
	ax4.text(0.6, 0.35, pipe_description[1], fontsize=10, color='blue')

	ax4.set_axis_off()

	for i in data:
		ax1.plot(*i)
		ax1.plot(*i, 'bo')
		for a, b in zip(*i):
			a, b = int(a), int(b)
			ax1.text(a, b, (a, b), ha='center', va='bottom', fontsize=8)
	
	ax1.text(data[0][0][0], data[0][1][0], shear_forces_description[0], ha='center', va='top', fontsize=8, color='red')
	ax1.text(data[1][0][0], data[1][1][0], shear_forces_description[1], ha='center', va='top', fontsize=8, color='red')

	gx = [i for i in range(rd[4])]
	ax5 = plt.subplot(plt.subplot(gs[2, :]))
	ax5.plot(gx, rd[3])
	plt.ylabel('Volume')
	plt.xlabel('Generation')

plt.show()
