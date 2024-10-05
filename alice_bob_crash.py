import numpy as np
import matplotlib.pyplot as plt

def rotate_txy(theta):
	"""
	Create a rotation matrix for (t,x,y) coordinates
	"""
	C = np.cos(theta)
	S = np.sin(theta)
	return np.array([[1,0,0],[0,C,S],[0,-S,C]])

def boostx_txy(beta):
	"""
	Lorentz boost matrix for the x direction in (t,x,y) coordinates.
	
	Meausre time in metres/c (units of 3.3ns)
	"""
	gamma = 1/np.sqrt(1-beta**2)
	return np.array([[gamma, -gamma*beta, 0],[-gamma*beta,gamma,0],[0,0,1]])
	
def boost_rot(theta,beta):
	"""
	Check of pen and paper math!
	"""
	C = np.cos(theta)
	S = np.sin(theta)
	gamma = 1/np.sqrt(1-beta**2)
	return np.array([[gamma,-gamma*beta*C,-gamma*beta*S],\
		[-gamma*beta*C,C**2*gamma+S**2, C*S*(gamma-1)],\
		[-gamma*beta*S, C*S*(gamma-1), S**2*gamma + C**2]])

#Define our time array, and our coordinates in Alice's frame
nt=4
t = np.linspace(-15,15,nt)
z = np.zeros(nt)
one = np.ones(nt)
C = np.cos(np.radians(10))
S = np.sin(np.radians(10))
A1 = np.array([t,z,z])
A2 = np.array([t,10*one,z])
B1 = np.array([t,0.5-0.9*C*t, 0.9*S*t])
B2 = np.array([t,9.5-0.9*C*t, 0.9*S*t])

#Now transform to Bob's frame
lorentz = boost_rot(np.radians(-10), -0.9)
A1p = np.dot(lorentz, A1)
A2p = np.dot(lorentz, A2)
B1p = np.dot(lorentz, B1)
B2p = np.dot(lorentz, B2)

#Plot the problem in Alice's frame.
plt.figure(1)
plt.clf()
plt.title('Alice Frame')
plt.plot(B1[2],B1[1], 'b')
plt.plot(B2[2],B2[1], 'b')
plt.plot(B1[2],B1[1], 'bo')
plt.plot(B2[2],B2[1], 'bo')
plt.plot(A1[2],A1[1], 'g*')
plt.plot(A2[2],A2[1], 'g*')
plt.xlabel('y(m)')
plt.ylabel('x(m)')
plt.tight_layout()

#Just plot the transformed coordinates! It is only 
#misleading because the positions of the line endpoints for Alice are
#actually at different times.
plt.figure(2)
plt.clf()
plt.title('Bob Frame??')
plt.plot(B1p[2],B1p[1], 'bo')
plt.plot(B2p[2],B2p[1], 'bo')
plt.plot(A1p[2],A1p[1], 'g:')
plt.plot(A2p[2],A2p[1], 'g:')
plt.xlabel('y(m)')
plt.ylabel('x(m)')
plt.tight_layout()

#Now lets interpolate Alice's spacecraft positions onto a fixed grid of times
#in Bob's frame. Interpolation is more complex than needed in this case, but
#allows us to consider accelerating frames in a similar way.
tint = np.linspace(-10,30,nt)
A1pi = np.empty_like(A1p)
A2pi = np.empty_like(A2p)
A1pi[0] = tint
A1pi[1] = np.interp(tint, A1p[0], A1p[1])
A1pi[2] = np.interp(tint, A1p[0], A1p[2])
A2pi[0] = tint
A2pi[1] = np.interp(tint, A2p[0], A2p[1])
A2pi[2] = np.interp(tint, A2p[0], A2p[2])

#Now plot Bob's frame, where we plot Alice's spacecraft at the same times 
#(relevant to Bob's instantaneous perception of the gap between spacecraft)
plt.figure(3)
plt.clf()
plt.title('Bob Frame - Alice Simultaneous')
plt.plot(B1p[2],B1p[1], 'bo')
plt.plot(B2p[2],B2p[1], 'bo')
plt.plot(A1pi[2],A1pi[1], 'g')
plt.plot(A2pi[2],A2pi[1], 'g')
plt.plot(A1pi[2],A1pi[1], 'g*')
plt.plot(A2pi[2],A2pi[1], 'g*')
plt.xlabel('y(m)')
plt.ylabel('x(m)')
plt.tight_layout()

