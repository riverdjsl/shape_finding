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
	dx = P2[0]-P1[0]
	dy = P2[1]-P1[1]
	l = np.sqrt(dx**2+dy**2)
	try:
		if dx > 0 and dy >= 0:
			return l, np.arctan(dy/dx)
		elif dy > 0 and dx <= 0:
			return l, np.arctan(abs(dx)/dy)+np.pi/2
		elif dx < 0 and dy <= 0:
			return l, np.arctan(dy/dx)+np.pi
		elif dy < 0 and dx >=0:
			return l, np.arctan(dx/abs(dy))+1.5*np.pi
		elif dy == 0 and dx == 0:
			return 0, 0
	except ZeroDivisionError:
		return 0, 0

def mesh(P1, P2, n=1):
	Xs = min(P1[0], P2[0])
	Xe = max(P1[0], P2[0])
	Ys = min(P1[1], P2[1])
	Ye = max(P1[1], P2[1])
	X = np.arange(Xs, Xe, n)
	Y = np.arange(Ys, Ye, n)
	X, Y = np.meshgrid(X, Y)
	sp = X.shape
	result = []
	for i in range(sp[0]):
		for j in range(sp[1]):
			result.append((X[i][j], Y[i][j]))
	return result

class Steel345():
	def __init__(s, t):
		s.t = t

	def fy(s):
		if s.t <= 16:
			return 305
		elif 16 < s.t <= 40:
			return 295
		else:
			return None

	def E(s):
		return 206e3

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

def fi(lambdai, f, E):
	alpha1, alpha2, alpha3 = 0.41, 0.986, 0.152
	lambdan = lambdai*np.sqrt(f/E)/np.pi
	if 0 < lambdan <= 0.215:
		return 1-alpha1*lambdan**2
	elif lambdan > 0.215:
		return ((alpha2+alpha3*lambdan+lambdan**2)-\
		np.sqrt((alpha2+alpha3*lambdan+lambdan**2)**2-4*lambdan**2))\
		/(2*lambdan**2)

def check(ratio, N, M, l, angle, u3, u2, E, f, W, I, A, i, r, t):
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

def NM(F, M0, P1, P2):
	gm = geometry(P1, P2)
	N = F*np.sin(gm[1])
	M = F*np.cos(gm[1])*gm[0]+M0
	return N, M, gm[0], gm[1],

def strainenergy(N, M, l, angle, u3, u2, E, f, W, I, A, i, r, t):
	i = np.sin(angle)
	j = np.cos(angle)
	part1 = (N*i)**2*l/(2*E*A)
	part2 = (N*j)**2*(l*r)**3*t*np.pi/(6*E*I**2)
	part3 = M**2*r**3*t*l*np.pi/(2*E*I**2)
	part4 = N*j*M*l**2*r**3*t*np.pi/(2*E*I**2)
	return part1+part2+part3+part4

def vars(pp, mt, mu):
	return {'u3': mu[0], 'u2': mu[1], 'E': mt.E(), 'f': mt.fy(), \
			'W': pp.W(), 'I': pp.I(), 'A': pp.A(), 'i': pp.i(), \
			'r': pp.r, 't': pp.t}

def pointlist4(x0, x1, x2, x3, x4, y):
	return [(x1, 0), (x4, 0), (x4, y), (x1, y), \
		(x1, y), (x2, y), (x3, y), (x4, y), (x0, 0)]

def qualify1(prec, P0, P1, P2, P3, A, B, C, D, PC):
	'''Mid branches moves arbitrarily.'''
	mesh1 = mesh(P0, P2, n=prec[0])
	for i in mesh1:
		if i[1] > 0:
			for j in mesh1:
				for k in mesh1:
					if A[1] > j[1] > i[1] and C[1] > k[1] > i[1]:
						if A[0] < j[0] <= i[0] <= k[0] < D[0]:
							yield [A, j], [B, j], [C, k], [D, k], \
								[j, i], [k, i], [i, PC]

def qualify2(prec, P0, P1, P2, P3, A, B, C, D, PC):
	'''The trunk top moves vertically.'''
	mesh1 = mesh(P0, P2, n=prec[0])
	for i in mesh1:
		if i[1] > 0 and i[0] == PC[0]:
			for j in mesh1:
				for k in mesh1:
					if A[1] > j[1] > i[1] and C[1] > k[1] > i[1]:
						if A[0] < j[0] <= i[0] <= k[0] < D[0]:
							yield [A, j], [B, j], [C, k], [D, k], \
								[j, i], [k, i], [i, PC]

def qualify3(prec, P0, P1, P2, P3, A, B, C, D, PC):
	'''The trunk top and the top of mid branch 1 move vertically.'''
	mesh1 = mesh(P0, P2, n=prec[0])
	for i in mesh1:
		if i[1] > 0 and i[0] == PC[0]:
			for j in mesh1:
				for k in mesh1:
					if A[1] > j[1] > i[1] and C[1] > k[1] > i[1]:
						if A[0] < j[0] < B[0] and \
							i[0] <=  k[0]  < D[0]:
							yield [A, j], [B, j], [C, k], [D, k], \
								[j, i], [k, i], [i, PC]

def qualify4(prec, P0, P1, P2, P3, A, B, C, D, PC):
	'''The trunk top and tops of mid branchs move vertically.'''
	mesh1 = mesh(P0, P2, n=prec[0])
	for i in mesh1:
		if i[1] > 0 and i[0] == PC[0]:
			for j in mesh1:
				for k in mesh1:
					if A[1] > j[1] > i[1] and C[1] > k[1] > i[1]:
						if A[0] < j[0] < B[0] and \
						C[0] < k[0] < D[0]:
							yield [A, j], [B, j], [C, k], [D, k], \
								[j, i], [k, i], [i, PC]

def qualify5(prec, P0, P1, P2, P3, A, B, C, D, PC):
	'''The trunk top and the top the 1st mid branch in one vertical line.'''
	mesh1 = mesh(P0, P2, n=prec[0])
	for i in mesh1:
		if i[1] > 0 and i[0] == PC[0]:
			for j in mesh1:
				if j [0] == PC[0]:
					for k in mesh1:
						if A[1] > j[1] > i[1] and C[1] > k[1] > i[1]:
							if i[0] < k[0] < D[0]:
								yield [A, j], [B, j], [C, k], [D, k], \
									[j, i], [k, i], [i, PC]

def qualify6(prec, P0, P1, P2, P3, A, B, C, D, PC):
	'''The trunk top and the tops of mid branches in one vertical line.'''
	mesh1 = mesh(P0, P2, n=prec[0])
	for i in mesh1:
		if i[1] > 0 and i[0] == PC[0]:
			for j in mesh1:
				if j [0] == PC[0]:
					for k in mesh1:
						if k[0] == PC[0]:
							if A[1] > j[1] > i[1] and C[1] > k[1] > i[1]:
									yield [A, j], [B, j], [C, k], [D, k], \
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

f = [qualify1, qualify2, qualify3, qualify4, qualify5, qualify6]
qualify = funcchoose(f)

datain = []
with open('input.txt') as f:
	for i in f.readlines():
		datain.append(trans_comma_str(i))

dataindic = {}
for j in datain:
	dicj = {}
	dicj[j[0]] = [float(k) for k in j[1:]]
	dataindic.update(dicj)

material = [Steel345(dataindic['topbranch'][1]), \
			Steel345(dataindic['middlebranch'][1]), \
			Steel345(dataindic['trunk'][1])]
pipes = [Pipe(*dataindic['topbranch']), \
		Pipe(*dataindic['middlebranch']), \
		Pipe(*dataindic['trunk'])]
mus = [dataindic['mutb'], dataindic['mumb'], dataindic['mut']]
vars1, vars2, vars3 = (vars(i, j, k) for i, j, k in zip(pipes, material, mus))
plist = pointlist4(*dataindic['geometry'])
F = dataindic['forces']
ratio = dataindic['ratio']
precision = dataindic['precision']

energies = []
collection = {}
for series_No, candidate in enumerate(qualify(precision, *plist)):
	collection[series_No] = candidate
	print('verifying: possible config '+str(series_No))

	energy = []
	checks = []

	dic1 = {}
	for j, k in zip(range(4), F):
		a = NM(k, 0, *candidate[j])
		dic1[j] = a

	dic2 = {}
	dic2[4] = NM(F[0]+F[1], dic1[0][1]+dic1[1][1], *candidate[4])
	dic2[5] = NM(F[2]+F[3], dic1[2][1]+dic1[3][1], *candidate[5])

	dic3 = {}
	dic3[6] = NM(sum(F), dic2[4][1]+dic2[5][1], *candidate[6])

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

final = collection[energies.index(min(energies))]

data = []
for i in final:
	xs = []
	ys = []
	for j in i:
		xs.append(j[0])
		ys.append(j[1])
	data.append([xs, ys])

for i in data:
	plt.plot(*i)
for i in data:
	plt.plot(*i, 'bo')
	for a, b in zip(*i):
		plt.text(a, b, (a, b), ha='center', va='bottom', fontsize=10)
plt.show()


