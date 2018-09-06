import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

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

'''
def NM(Nv, M0, lv):
	l = np.sqrt(lv.dot(lv))
	elv = lv/l
	N = np.sqrt(Nv.dot(Nv))
	eNv = Nv/N
	Mv = np.cross(lv, Nv)
	eMv = Mv/np.sqrt(Mv.dot(Mv))
	eMv = np.nan_to_num(eMv)

	try:
		thitar = np.arccos(lv.dot(Nv)/(l*N))
	except Exception:
		thitar = 0
	P = N*np.cos(thitar)*elv
	M = N*np.sin(thitar)*np.sqrt(lv.dot(lv))*eMv+M0
	return P, M, l, thitar
'''

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


def check(ratio, N, l, u3, u2, E, f, W, I, A, i, r, t):

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
	
	if condition1 <= ratio and condition2 <= ratio:
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


def chromosome_translator(individual, tmax, A, B, PC):

	def coord(chromosome, a, b):
		return a+(b-a)*DEC(chromosome)/511

	yi = coord(individual[0], 0, A[1])
	i = (PC[0], yi)

	xj = coord(individual[1], A[0], PC[0])
	yj = coord(individual[2], i[1], A[1])
	j = (xj, yj)

	t1 = coord(individual[3], 5, tmax[0])
	t2 = 2*t1

	return ([A, j], [B, j], [j, i], [i, PC]), (t1, t2), (i, j)


#sfitness points: the minimun value gets the highest point.	
def point1(rou, v, TF):
	'''Fitness based on total volumes!'''
	if TF is True:
		return 5-v/1e9
	else:
		return 0


def point2(rou, v, TF):
	'''Fitness based on energy density.'''
	if TF is True:
		return rou
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

	p_num = int(dataindic['population'][0])
	g_num = int(dataindic['generation'][0])

	F_direction = np.array([0, -1, 0])

	def energy(individual, F):

		individual_values = chromosome_translator(individual, tmax, *plist)
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

		dic1 = {}
		for j, k in zip(range(8), F):
			a = NM(k*F_direction, np.zeros(3), np.array(candidate[j][0])-np.array(candidate[j][1]))
			dic1[j] = a

		dic2 = {}
		dic2[8] = NM((F[0]+F[1])*F_direction, dic1[0][1]+dic1[1][1], np.array(candidate[8][0])-np.array(candidate[8][1]))
		dic2[10] = NM((F[4]+F[5])*F_direction, dic1[4][1]+dic1[5][1], np.array(candidate[10][0])-np.array(candidate[10][1]))
		dic2[9] = NM((F[2]+F[3])*F_direction, dic1[2][1]+dic1[3][1], np.array(candidate[9][0])-np.array(candidate[9][1]))
		dic2[11] = NM((F[6]+F[7])*F_direction, dic1[6][1]+dic1[7][1], np.array(candidate[11][0])-np.array(candidate[11][1]))

		dic3 = {}
		dic3[12] = NM(sum(F)*F_direction, dic2[8][1]+dic2[9][1]+dic2[10][1]+dic2[11][1], np.array(candidate[12][0])-np.array(candidate[12][1]))

		for i, j in dic1.items():
			energyi.append(strainenergy(*j, **vars1)[0])
			volumni.append(strainenergy(*j, **vars1)[1])
			checks.append(check(ratio[0], *j, **vars1))
		for i, j in dic2.items():
			energyi.append(strainenergy(*j, **vars2)[0])
			volumni.append(strainenergy(*j, **vars1)[1])
			checks.append(check(ratio[1], *j, **vars2))
		for i, j in dic3.items():
			energyi.append(strainenergy(*j, **vars3)[0])
			volumni.append(strainenergy(*j, **vars1)[1])
			checks.append(check(ratio[2], *j, **vars3))

		energy = sum(energyi)
		volumn = sum(volumni)
		density = energy/volumn
		
		return density, volumn, all(checks)

	population = []
	for i in range(p_num):
		population.append(individual_generator(np.random.rand(16)))

	generation = 0
#	energy_g = list(map(energy, population, [forces]*len(population)))

	animate = []
	volumnmins = []
	while generation <= g_num:

		print('\n\n\nGeneration {}.\n'.format(generation))

		volumn_g = list(map(energy, population, [forces]*len(population)))
		pts = [point(i, j, k) for i, j, k in volumn_g]

		if generation > 0.6*g_num and abs(np.average(pts)-max(pts)) <= 0.005:
			print("Ending the run 'cause convergence almost achieved!")
			break

		print('\n\n'.join(list('  '.join(list('{:>5.2f}'.format(i) for i in pts)[j: j+5]) for j in range(0, len(pts), 5))))

		pama = select(population, pts)
		kids = reproduction(pama)
		v_kids = list(map(energy, kids, [forces]*len(kids)))
		pts_k = [point(i, j, k) for i, j, k in v_kids]

		for i in sorted(pts[:])[:4]:
			for j in pts_k[:]:
				if j >= i:
					population[pts.index(i)] = kids[pts_k.index(j)]
					pts_k.pop(pts_k.index(j))

		if abs(np.average(pts)-max(pts)) > 1:
			population = population_mutation(population)
		else:
			population = population_mutation(population, possibility=0.05)

		volumnmins.append(min(volumn_g)[1])
		print('\nCurrent Volumn: {:.1f}.'.format(min(volumn_g)[1]))
		individual_em = population[volumn_g.index(min(volumn_g))]

		animate_i = [[]]*5
		if generation%1 == 0:
			stepdata = chromosome_translator(individual_em, tmax, *plist)
			animate.append(stepdata[2])

		generation += 1

	return chromosome_translator(individual_em, tmax, *plist), dataindic['D_over_t'][0], min(volumn_g)[1], volumnmins, generation, ratio, animate




rd = rawdata('input_GA.txt')

data = []
for i in rd[0][0]:
	xs = []
	ys = []
	zs = []
	for j in i:
		xs.append(j[0])
		ys.append(j[1])
		zs.append(j[2])
	data.append([xs, ys, zs])

pipe_description = []
for i,j in zip(rd[0][1], rd[5]):
	pipe_description.append('D = {:.1f}, t = {:.1f}, σ/f ≈ {:.1f}'.format(i*rd[1], i, j))

fig1 = plt.figure('fig1')
ax1 = plt.subplot(321)
ax1.axis('equal')
ax1.set_title('Overlook View')
plt.ylabel('Z')
ax2 = plt.subplot(323)
ax2.axis('equal')
ax2.set_title('Front View')
plt.ylabel('Y')
plt.xlabel('X')
ax3 = plt.subplot(324)
ax3.axis('equal')
ax3.set_title('Side View')
plt.xlabel('Z')
ax4 = plt.subplot(322)
ax4.text(0.0, 0.6, 'The total volume:', fontsize=10)
ax4.text(0.1, 0.5, int(rd[2]), fontsize=10, color='red')
ax4.text(0.4, 0.6, 'The uper branches:', fontsize=10)
ax4.text(0.4, 0.5, pipe_description[0], fontsize=10, color='blue')
ax4.text(0.4, 0.4, 'The middle branches:', fontsize=10)
ax4.text(0.4, 0.3, pipe_description[1], fontsize=10, color='blue')
ax4.text(0.4, 0.2, 'The trunk:', fontsize=10)
ax4.text(0.4, 0.1, pipe_description[2], fontsize=10, color='blue')
ax4.set_axis_off()

for i in data:
	ax1.plot(*i[0:3:2])
	ax2.plot(*i[:2])
	ax3.plot(*i[::-1][:2])

	ax1.plot(*i[0:3:2], 'bo')
	for a, b in zip(*i[0:3:2]):
		a, b = int(a), int(b)
		ax1.text(a, b, (a, b), ha='center', va='bottom', fontsize=10)
	ax2.plot(*i[:2], 'bo')
	for c, d in zip(*i[:2]):
		c, d = int(c), int(d)
		ax2.text(c, d, (c, d), ha='center', va='bottom', fontsize=10)
	ax3.plot(*i[::-1][:2], 'bo')
	for e, f in zip(*i[::-1][:2]):
		e, f = int(e), int(f)
		ax3.text(e, f, (e, f), ha='center', va='bottom', fontsize=10)

gx = [i for i in range(rd[4])]
ax5 = plt.subplot(313)
ax5.plot(gx, rd[3])
plt.ylabel('Volume')
plt.xlabel('Generation')


joint_num = len(rd[6][0])
joint_dic = {}
for i in range(joint_num):
	joint_dic[str(i)] = []

for j in rd[6]:
	for i in range(joint_num):
		joint_dic[str(i)].append(j[i])

joint_dic2 = {}
for i, j in joint_dic.items():
	xistep = []
	yistep = []
	zistep = []
	for k in j:
		xistep.append(k[0])
		yistep.append(k[2])
		zistep.append(k[1])
	joint_dic2[i] = np.array([xistep, yistep, zistep])

fig2 = plt.figure('Animation')
ax = Axes3D(fig2)


def update_lines(num, dataLines, lines):
	for line, data in zip(lines, dataLines):
		line.set_data(data[0:2, num-50:num])
		line.set_3d_properties(data[2, num-50:num])
	return lines

marks = ['r^', 'bx', 'g*', 'yD', 'o']

data = list(joint_dic2.values())
lines = [ax.plot(dat[0], dat[1], dat[2], marker)[0] for dat, marker in zip(data, marks)]

line_ani = FuncAnimation(fig2, update_lines, len(data[0][0]), fargs=(data, lines), interval=1, blit=True)

plt.show()
