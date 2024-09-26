import numpy as np
import control as cnt
import massParam as P

class ctrlStateFeedbackIntegrator:
    def __init__(self):
        #  tuning parameters
        tr = 2.5
        zeta = 0.95
        integrator_pole = -1.2

        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        A = np.array([[0.0, 1.0],
                      [-P.k / P.m, -P.b / P.m]])
        B = np.array([[0.0],
                      [1.0 / P.m]])
        C = np.array([[1.0, 0.0]])

        # form augmented system
        A1 = np.array([[0.0, 1.0, 0.0],
                       [-P.k / P.m, -P.b / P.m, 0.0],
                       [-1.0, 0.0, 0.0]])
        B1 = np.array([[0.0],
                       [1.0 / P.m],
                       [0.0]])
        
        # gain calculation
        wn = 0.5*np.pi/(tr*np.sqrt(1-zeta**2)) # natural frequency for when zeta is not equal to 0.707
        des_char_poly = np.convolve([1, 2*zeta*wn, wn**2],
                                    [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 3:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)  #the cnt.acker function also works
            self.K = K1[0][0:2]
            self.ki = K1[0][2]
        print('K: ', self.K)
        print('ki: ', self.ki)

        #--------------------------------------------------
        # variables to implement integrator
        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample

    def update(self, z_r, x):
        z = x[0][0]
        z_dot = x[1][0]

        # integrate error
        error = z_r - z
        self.integrator = self.integrator \
                          + (P.Ts / 2.0) * (error + self.error_d1)
        self.error_d1 = error

        # Compute the state feedback controller
        force_unsat = -self.K @ x - self.ki*self.integrator

        # compute total torque
        force = saturate(force_unsat[0], P.F_max)
        
        # integrator anti-windup
        if self.ki != 0.0:
            self.integrator = self.integrator + P.Ts/self.ki*(force-force_unsat)
        return force


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u

