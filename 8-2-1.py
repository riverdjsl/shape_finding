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
	

def mesh(P1, P2, n=1):
	Xs = min(P1[0], P2[0])
	Xe = max(P1[0], P2[0])
	Ys = min(P1[1], P2[1])
	Ye = max(P1[1], P2[1])
	Zs = min(P1[2], P2[2])
	Ze = max(P1[2], P2[2])	
	X = np.arange(Xs, Xe, n)
	Y = np.arange(Ys, Ye, n)
	Z = np.arange(Zs, Ze, n)
	X, Y = np.meshgrid(X, Y)
	sp = X.shape
	result = []
	for k in Z:
		for i in range(sp[0]):
			for j in range(sp[1]):
				result.append((X[i][j], Y[i][j], k))
	return result


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


def mtgenerater(t):
	return [Steel345(t), Steel345(2*t),\
			Steel345(2*np.sqrt(2)*t)]


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


def pipegenerater(t, factor):
	return [Pipe(t*factor, t), Pipe(2*t*factor, 2*t),\
			Pipe(2*np.sqrt(2)*t*factor, 2*np.sqrt(2)*t)]


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
	return [(x1a, 0, z1a), (x4, y4, z4), \
		(x1, y1, z1), (x2, y2, z2), (x3, y3, z3), (x4, y4, z4), \
		(x1a, y1a, z1a), (x2a, y2a, z2a), (x3a, y3a, z3a), \
		(x4a, y4a, z4a), (x0, 0, 0)]


def qualify1(prec, P0, P1, A, B, C, D, Aa, Ba, Ca, Da, PC):
	'''Mid branches moves arbitrarily.'''
	mesh1 = mesh(P0, P1, n=prec[0])
	for i in mesh1:
		if i[1] > 0:
			for j in mesh1:
				for k in mesh1:
					if A[1] > j[1] > i[1] and C[1] > k[1] > i[1]:
						if A[0] < j[0] < i[0] < k[0] < D[0]:
							yield [A, j], [B, j], [C, k], [D, k], \
								[Aa, j], [Ba, j], [Ca, k], [Da, k], \
								[j, i], [k, i], [i, PC]
			

def qualify2(prec, P0, P1, A, B, C, D, Aa, Ba, Ca, Da, PC):
	'''The trunk top moves vertically.'''
	mesh1 = mesh(P0, P1, n=prec[0])
	for i in mesh1:
		if i[1] > 0 and i[0] == PC[0] and i[2] == 0:
			for j in mesh1:
				for k in mesh1:
					if A[1] > j[1] > i[1] and C[1] > k[1] > i[1]:
						if A[0] < j[0] < i[0] < k[0] < D[0]:
							yield [A, j], [B, j], [C, k], [D, k], \
								[Aa, j], [Ba, j], [Ca, k], [Da, k], \
								[j, i], [k, i], [i, PC]
								

def funcchoose(funclist):
	for i, j in enumerate(funclist):
		print('Method ', i, ':', j.__doc__)
	while True:
		try:
			fn = input("Now choose one:")
			for i, j in enumerate(funclist):
				if i == int(fn):
					return j
		except Exception:
			continue

f = [qualify1, qualify2]
qualify = funcchoose(f)

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

	material = mtgenerater(dataindic['t'][0])
	pipes = pipegenerater(dataindic['t'][0], dataindic['D_over_t'][0])
	mus = [dataindic['mutb'], dataindic['mumb'], dataindic['mut']]
	vars1, vars2, vars3 = (vars(i, j, k) \
		for i, j, k in zip(pipes, material, mus))
	plist = pointlist8(*dataindic['geometry'])
	F = dataindic['forces']
	ratio = dataindic['ratio']
	precision = dataindic['precision']

	F_direction = np.array([0, -1, 0])
	energies = []
	collection = {}
	for series_No, candidate in enumerate(qualify(precision, *plist)):
		collection[series_No] = candidate
		print('verifying: possible config '+str(series_No))

		energy = []
		checks = []

		dic1 = {}
		for j, k in zip(range(8), F):
			a = NM(k*F_direction, np.zeros(3), np.array(candidate[j][0])-np.array(candidate[j][1]))
			dic1[j] = a

		dic2 = {}
		dic2[8] = NM((F[0]+F[1]+F[4]+F[5])*F_direction, dic1[0][1]+dic1[1][1]+dic1[4][1]+dic1[5][1], np.array(candidate[8][0])-np.array(candidate[8][1]))
		dic2[9] = NM((F[2]+F[3]+F[6]+F[7])*F_direction, dic1[2][1]+dic1[3][1]+dic1[6][1]+dic1[7][1], np.array(candidate[9][0])-np.array(candidate[9][1]))
	
		dic3 = {}
		dic3[10] = NM(sum(F)*F_direction, dic2[8][1]+dic2[9][1], np.array(candidate[10][0])-np.array(candidate[10][1]))

		for i, j in dic1.items():
			energy.append(strainenergy(*j, **vars1))
			checks.append(check(ratio[0], *j, **vars1))
		for i, j in dic2.items():
			energy.append(strainenergy(*j, **vars2))
			checks.append(check(ratio[1], *j, **vars2))
		for i, j in dic3.items():
			energy.append(strainenergy(*j, **vars3))
			checks.append(check(ratio[2], *j, **vars3))

		if all(checks):
			energies.append(sum(energy))
	try:
		return collection[energies.index(min(energies))], pipes

	except Exception:
		print('error')

data = []
rd = rawdata('input8.txt')
for i in rd[0]:
	xs = []
	ys = []
	zs = []
	for j in i:
		xs.append(j[0])
		ys.append(j[1])
		zs.append(j[2])
	data.append([xs, ys, zs])
	
print(data)

for i in rd[1]:
	print('D = {:.1f}, t = {:.1f}'.format(i.D, i.t))



ax1 = plt.subplot(221)
ax1.axis('equal')
ax2 = plt.subplot(223)
ax2.axis('equal')
ax3 = plt.subplot(224)
ax3.axis('equal')
for i in data:
	ax1.plot(*i[0:3:2])
	ax2.plot(*i[:2])
	ax3.plot(*i[::-1][:2])

	ax1.plot(*i[0:3:2], 'bo')
	for a, b in zip(*i[0:3:2]):
		ax1.text(a, b, (a, b), ha='center', va='bottom', fontsize=10)
	ax2.plot(*i[:2], 'bo')
	for c, d in zip(*i[:2]):
		ax2.text(c, d, (c, d), ha='center', va='bottom', fontsize=10)
	ax3.plot(*i[::-1][:2], 'bo')
	for e, f in zip(*i[::-1][:2]):
		ax3.text(e, f, (e, f), ha='center', va='bottom', fontsize=10)

plt.show()


