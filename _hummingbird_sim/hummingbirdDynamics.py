import numpy as np 
import hummingbirdParam as P




class HummingbirdDynamics:
   def __init__(self, alpha=0.0):
       # Initial state conditions
       self.state = np.array([
           [P.phi0],  # roll angle
           [P.theta0],  # pitch angle
           [P.psi0],  # yaw angle
           [P.phidot0],  # roll rate
           [P.thetadot0],  # pitch rate
           [P.psidot0],  # yaw rate
       ])

       # vary the actual physical parameters
       self.ell1 = P.ell1 * (1.+alpha*(2.*np.random.rand()-1.))
       self.ell2 = P.ell2 * (1.+alpha*(2.*np.random.rand()-1.))
       self.ell3x = P.ell3x * (1.+alpha*(2.*np.random.rand()-1.))
       self.ell3y = P.ell3y * (1.+alpha*(2.*np.random.rand()-1.))
       self.ell3z = P.ell3z * (1.+alpha*(2.*np.random.rand()-1.))
       self.ellT = P.ellT * (1.+alpha*(2.*np.random.rand()-1.))
       self.d = P.d * (1.+alpha*(2.*np.random.rand()-1.))
       self.m1 = P.m1 * (1.+alpha*(2.*np.random.rand()-1.))
       self.m2 = P.m2 * (1.+alpha*(2.*np.random.rand()-1.))
       self.m3 = P.m3 * (1.+alpha*(2.*np.random.rand()-1.))
       self.J1x = P.J1x * (1.+alpha*(2.*np.random.rand()-1.))
       self.J1y = P.J1y * (1. + alpha * (2. * np.random.rand() - 1.))
       self.J1z = P.J1z * (1. + alpha * (2. * np.random.rand() - 1.))
       self.J2x = P.J2x * (1.+alpha*(2.*np.random.rand()-1.))
       self.J2y = P.J2y * (1. + alpha * (2. * np.random.rand() - 1.))
       self.J2z = P.J2z * (1. + alpha * (2. * np.random.rand() - 1.))
       self.J3x = P.J3x * (1.+alpha*(2.*np.random.rand()-1.))
       self.J3y = P.J3y * (1. + alpha * (2. * np.random.rand() - 1.))
       self.J3z = P.J3z * (1. + alpha * (2. * np.random.rand() - 1.))
       self.km = P.km * (1. + alpha * (2. * np.random.rand() - 1.))
   def update(self, u: np.ndarray):
       # This is the external method that takes the input u at time
       # t and returns the output y at time t.
       # saturate the input force
       u = saturate(u, P.torque_max)
       self.rk4_step(u)  # propagate the state by one time sample
       y = self.h()  # return the corresponding output
       return y

   def f(self, state: np.ndarray, pwms: np.ndarray):
       # Return xdot = f(x,u)
       phidot = state[3][0]
       thetadot = state[4][0]
       psidot = state[5][0]
       pwm_left = pwms[0][0]
       pwm_right = pwms[1][0]

       # The equations of motion go here
       M = self._M(state)
       C = self._C(state)
       partialP = self._partialP(state)

       force = self.km * (pwm_left + pwm_right)
       torque = self.d * self.km * (pwm_left - pwm_right)
       tau = self._tau(state, force, torque)
       B = self._B()

       qddot = np.linalg.inv(M) @ (-C - partialP + tau - B @ state[3:6])
       
       phiddot = qddot[0][0]
       thetaddot = qddot[1][0]
       psiddot = qddot[2][0]
       
       # build xdot and return
       xdot = np.array([[phidot],
                        [thetadot],
                        [psidot],
                        [phiddot],
                        [thetaddot],
                        [psiddot]])
       return xdot

   def h(self):
       # FIXME Fill in this function
       # return y = h(x)
       phi = self.state[0][0]
       theta = self.state[1][0]
       psi = self.state[2][0]
       y = np.array([[phi], [theta], [psi]])
       return y

   def rk4_step(self, u: np.ndarray):
       # Integrate ODE using Runge-Kutta RK4 algorithm
       F1 = self.f(self.state, u)
       F2 = self.f(self.state + P.Ts / 2 * F1, u)
       F3 = self.f(self.state + P.Ts / 2 * F2, u)
       F4 = self.f(self.state + P.Ts * F3, u)
       self.state = self.state + P.Ts / 6 * (F1 + 2*F2 + 2*F3 + F4)

   def _M(self, state: np.ndarray):
       # FIXME Fill in this function
       phi = state[0][0]
       theta = state[1][0]
       psi = state[2][0]
       phidot = state[3][0]
       thetadot = state[4][0]
       psidot = state[5][0]

       ctheta = np.cos(theta)
       cphi = np.cos(phi)
       cpsi = np.cos(psi)

       stheta = np.sin(theta)
       sphi = np.sin(phi)
       spsi = np.sin(psi)



       # Fill out M22, M23, and M33
       M22 = self.m1*self.ell1**2+ \
           self.m2*self.ell2**2 \
               +self.J2y+ \
                   self.J1y*cphi**2 \
                   +self.J1z*sphi**2
       M23 = (self.J1y-self.J1z)*sphi*cphi*ctheta
       M33 = (self.m1*self.ell1**2+self.m2*self.ell2**2+self.J2z+self.J1y*sphi**2+self.J1z*cphi**2)*ctheta**2+(self.J1x+self.J2x)*stheta**2+self.m3*(self.ell3x**2+self.ell3y**2)+self.J3z

       # Return the M matrix
       return np.array([[self.J1x, 0, -self.J1x*stheta],
                     [0, M22, M23],
                     [-self.J1x*stheta, M23, M33]
                     ])

   def _C(self, state: np.ndarray):
       # FIXME Fill in this function
       #extact any necessary variables from the state
       phi = state[0][0]
       theta = state[1][0]
       psi = state[2][0]
       phidot = state[3][0]
       thetadot = state[4][0]
       psidot = state[5][0]

       ctheta = np.cos(theta)
       cphi = np.cos(phi)
       cpsi = np.cos(psi)

       stheta = np.sin(theta)
       sphi = np.sin(phi)
       spsi = np.sin(psi)

       N33 = 2*(self.J1x+self.J2x\
                -self.m1*self.ell1**2\
                   -self.m2*self.ell2**2\
                       -self.J2z\
                -self.J1y*sphi**2\
                   -self.J1z*cphi**2)\
                       *stheta*ctheta
       
       C1 = (self.J1y-self.J1z)*sphi*cphi*\
                      (thetadot**2-ctheta**2*psidot**2)\
                      +((self.J1y-self.J1z)*(cphi**2-sphi**2)-self.J1x)\
                       *ctheta*thetadot*psidot

       C2 = 2*(self.J1z-self.J1y)*sphi*cphi*phidot*thetadot+\
                      ((self.J1y-self.J1z)*(cphi**2-sphi**2)+self.J1x)\
                       *ctheta*phidot*psidot-0.5*N33*psidot**2

       C3 = thetadot**2*(self.J1z-self.J1y)*sphi*cphi*stheta\
                       +((self.J1y-self.J1z)*(cphi**2-sphi**2)\
                         -self.J1x)*ctheta*phidot*thetadot\
                           +(self.J1z-self.J1y)*sphi*cphi*stheta*thetadot**2\
                               +2*(self.J1y-self.J1z)*sphi*cphi*phidot*psidot\
                                   +2*(-self.m1*self.ell1**2-self.m2*self.ell2**2-self.J2z\
                                       +self.J1x\
                                           +self.J2x+\
                                               self.J1y*sphi**2\
                                                   +self.J1z*sphi**2)*stheta\
                                                       *ctheta*thetadot*psidot

       # Return the C matrix
       return np.array([[C1],
                     [C2],
                     [C3],
                    ])
       
   def _partialP(self, state: np.ndarray):
       # FIXME Fill in this function
       #extact any necessary variables from the state
       phi = state[0][0]
       theta = state[1][0]
       psi = state[2][0]
       phidot = state[3][0]
       thetadot = state[4][0]
       psidot = state[5][0]

       ctheta = np.cos(theta)
       cphi = np.cos(phi)
       cpsi = np.cos(psi)

       stheta = np.sin(theta)
       sphi = np.sin(phi)
       spsi = np.sin(psi)

       # Return the partialP array
       return np.array([[0],
                            [(self.m1*self.ell1+self.m2*self.ell2)*P.g*ctheta],
                            [0],
                           ])
   
   def _tau(self, state: np.ndarray, force: float, torque: float):
       """
       Returns the tau matrix as defined in the hummingbird manual.

       Parameters
       ----------
       state : numpy.ndarray
           The state of the hummingbird. Contains phi, theta, psi, and their derivatives.
       force : float
           force = (fl + fr). e.g. the second element of the tau matrix becomes
           lT * force * cos(phi) using the above definition.
       torque : float
           torque = d(fl - fr). e.g. the first element of teh tau matrix just
           becomes torque, using the definition above.

       """
       # FIXME Fill in this function
       #extract any necessary variables from the state

       

       phi = state[0][0]
       theta = state[1][0]
       psi = state[2][0]
       phidot = state[3][0]
       thetadot = state[4][0]
       psidot = state[5][0]

       ctheta = np.cos(theta)
       cphi = np.cos(phi)
       cpsi = np.cos(psi)

       stheta = np.sin(theta)
       sphi = np.sin(phi)
       spsi = np.sin(psi)

       # Return the tau matrix
       return np.array([[torque],
                       [self.ellT*force*cphi],
                       [self.ellT*force*ctheta*sphi-torque*stheta]])
   def _B(self):
       # FIXME Fill in this function
       # This needs no variables from the state
       B = 0.001 * np.identity(3, dtype = float)
       
       # Return the B matrix
       return B


def saturate(u: np.ndarray, limit: float):
   for i in range(0, u.shape[0]):
       if abs(u[i][0]) > limit:
           u[i][0] = limit * np.sign(u[i][0])
   return u
