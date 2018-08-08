import numpy as np
import matplotlib.pyplot as plt

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


def geometry(P1, P2):
	return np.array(P2)-np.array(P1)


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
	return [Steel345(ts[0]), Steel345(ts[1]), Steel345(ts[2])]
			

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
	return [Pipe(ts[0]*factor, ts[0]), Pipe(ts[1]*factor, ts[1]), \
			Pipe(ts[2]*factor, ts[2])]


def fi(lambdai, f, E):
	alpha1, alpha2, alpha3 = 0.41, 0.986, 0.152
	lambdan = lambdai*np.sqrt(f/E)/np.pi
	if 0 < lambdan <= 0.215:
		return 1-alpha1*lambdan**2
	elif lambdan > 0.215:
		return ((alpha2+alpha3*lambdan+lambdan**2)-\
		np.sqrt((alpha2+alpha3*lambdan+lambdan**2)**2-4*lambdan**2))\
		/(2*lambdan**2)


def check(ratio, Nv, Mv, lv, angle, u3, u2, E, f, W, I, A, i, r, t):
	N = np.sqrt(Nv.dot(Nv))
	M = np.sqrt(Mv.dot(Mv))
	l = lv
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
	part1 = N/(fix*A*f)
	part2 = 0.6*M/(1.15*W*(1-0.8*N/Nex)*f)
	part3 = 0.7*M/(W*f)
	part4 = N/(A*f)
	part5 = M/(1.15*W*f)
	condition1 = abs(part4+part5)
	condition2 = abs(part1+part2)
	condition3 = abs(part1+part3)
	if condition1 < ratio and condition2 < ratio and condition3 < ratio:
		return True
	else:
		return False


def strainenergy(Nv, Mv, lv, angle, u3, u2, E, f, W, I, A, i, r, t):
	N = np.sqrt(Nv.dot(Nv))
	M = np.sqrt(Mv.dot(Mv))
	l = lv
	i = np.cos(angle)
	j = np.sin(angle)
	part1 = (N*i)**2*l/(2*E*A)
	part2 = (N*j)**2*(l*r)**3*t*np.pi/(6*E*I**2)
	part3 = M**2*r**3*t*l*np.pi/(2*E*I**2)
	part4 = N*j*M*l**2*r**3*t*np.pi/(2*E*I**2)
	return part1+(part2+part3+part4)


def vars(pp, mt, mu):
	return {'u3': mu[0], 'u2': mu[1], 'E': mt.E(), 'f': mt.fy(), \
			'W': pp.W(), 'I': pp.I(), 'A': pp.A(), 'i': pp.i(), \
			'r': pp.r, 't': pp.t}


def pointlist8(x0, x1, x2, x3, x4, x1a, x2a, x3a, x4a, y1, y2, y3, y4, \
				y1a, y2a, y3a, y4a, z1, z2, z3, z4, z1a, z2a, z3a, z4a):
	return [(x1, y1, z1), (x2, y2, z2), (x3, y3, z3), (x4, y4, z4), \
		(x1a, y1a, z1a), (x2a, y2a, z2a), (x3a, y3a, z3a), \
		(x4a, y4a, z4a), (x0, 0, 0)]


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

	
def chromosome_translator(individual, tmax, A, B, C, D, Aa, Ba, Ca, Da, PC):
		
	def coord(chromosome, a, b):
		return a+(b-a)*DEC(chromosome)/511
		
	yi = coord(individual[0], 0, min(A[1], B[1], C[1], D[1], Aa[1], Ba[1], Ca[1], Da[1]))
	i = (PC[0], yi, 0)

	xj = coord(individual[1], A[0], PC[0])
	yj = coord(individual[2], i[1], min(A[1], B[1]))
	zj = coord(individual[3], 0, max(A[2], B[2]))
	j = (xj, yj, zj)

	xja = coord(individual[4], Aa[0], PC[0])
	yja = coord(individual[5], i[1], min(Aa[1], Ba[1]))
	zja = coord(individual[6], 0, min(Aa[2], Ba[2]))
	ja = (xja, yja, zja)
	
	xk = coord(individual[7], PC[0], D[0])
	yk = coord(individual[8], i[1], min(C[1], D[1]))
	zk = coord(individual[9], 0, max(C[2], D[2]))
	k = (xk, yk, zk)
	
	xka = coord(individual[10], PC[0], Da[0])
	yka = coord(individual[11], i[1], min(Ca[1], Da[1]))
	zka = coord(individual[12], 0, min(Ca[2], Da[2]))
	ka = (xka, yka, zka)

	t1 = coord(individual[13], 5, tmax[0])
	t2 = coord(individual[14], 5, tmax[1])
	t3 = coord(individual[15], 5, tmax[2])
	
	return ([A, j], [B, j], [C, k], [D, k], [Aa, ja], [Ba, ja], [Ca, ka], [Da, ka], [j, i], [k, i], [ja, i], [ka, i], [i, PC]), (t1, t2, t3)


#sfitness points: the minimun value gets the highest point.
def point(value, TF):
	if TF is True:
		return 19-9*value/5e7
	else:
		return 0


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
	if chance >= 0.5:
		ch = list(chromosome)
		for i, gene in enumerate(ch):
			if  i == mutation_point:
				ch[i] = str((1-int(gene)))
		return ''.join(ch)
	else:
		return chromosome


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


	mus = [dataindic['mutb'], dataindic['mumb'], dataindic['mut']]
	plist = pointlist8(*dataindic['geometry'])
	forces = dataindic['forces']
	ratio = dataindic['ratio']
	tmax = dataindic['tmax']
	
	p_num = int(dataindic['population'][0])
	g_num = int(dataindic['generation'][0])
	
	F_direction = np.array([0, -1, 0])

	def energy(individual, F):
		
		candidate = chromosome_translator(individual, tmax, *plist)[0]
		pipets = chromosome_translator(individual,  tmax, *plist)[1]
		
		material = mtgenerater(pipets)
		pipes = pipegenerater(pipets, dataindic['D_over_t'][0])
		
		vars1, vars2, vars3 = (vars(i, j, k) \
		for i, j, k in zip(pipes, material, mus))
		
		energyi = []
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
			energyi.append(strainenergy(*j, **vars1))
			checks.append(check(ratio[0], *j, **vars1))
		for i, j in dic2.items():
			energyi.append(strainenergy(*j, **vars2))
			checks.append(check(ratio[1], *j, **vars2))
		for i, j in dic3.items():
			energyi.append(strainenergy(*j, **vars3))
			checks.append(check(ratio[2], *j, **vars3))

		energy = sum(energyi)
		return energy, all(checks)

	population = []	
	for i in range(p_num):
		population.append(individual_generator(np.random.rand(16)))

	generation = 0
	energy_g = list(map(energy, population, [forces]*len(population)))
	
	energymins = []
	while generation <= g_num:
		print(generation)
		pts = [point(i, j) for i, j in energy_g]
		print(pts)
		pama = select(population, pts)
		kids = reproduction(pama)
		e_kids = list(map(energy, kids, [forces]*len(kids)))
		pts_k = [point(i, j) for i, j in e_kids]
		for i in sorted(pts[:])[:4]:
			for j in pts_k[:]:
				if j >= i:
					population[pts.index(i)] = kids[pts_k.index(j)]
					pts_k.pop(pts_k.index(j))
		energy_g = list(map(energy, population, [forces]*len(population)))
		energymins.append(min(energy_g)[0])
		individual_em = population[energy_g.index(min(energy_g))]
		generation += 1

	return chromosome_translator(individual_em, tmax, *plist), dataindic['D_over_t'][0], min(energy_g), energymins, generation, 


data = []
rd = rawdata('input_GA.txt')
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
for i in rd[0][1]:
	pipe_description.append('D = {:.1f}, t = {:.1f}'.format(i*rd[1], i))

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
ax4.text(0.0, 0.6, 'The total energy:', fontsize=10)
ax4.text(0.1, 0.5, int(rd[2][0]), fontsize=10, color='red')
ax4.text(0.4, 0.6, 'The uper branches:', fontsize=10)
ax4.text(0.5, 0.5, pipe_description[0], fontsize=10, color='blue')
ax4.text(0.4, 0.4, 'The middle branches:', fontsize=10)
ax4.text(0.5, 0.3, pipe_description[1], fontsize=10, color='blue')
ax4.text(0.4, 0.2, 'The trunk:', fontsize=10)
ax4.text(0.5, 0.1, pipe_description[2], fontsize=10, color='blue')
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
plt.ylabel('Total Energy')
plt.xlabel('Generation')

plt.show()



