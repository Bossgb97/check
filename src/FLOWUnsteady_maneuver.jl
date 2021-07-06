#=##############################################################################
# DESCRIPTION
    Types defining maneuvers of flight vehicles.

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Oct 2019
  * License   : MIT
=###############################################################################

################################################################################
# ABSTRACT MANEUVER TYPE
################################################################################
"""
    `AbstractManeuver{N, M}`

`N` indicates the number of tilting systems in this maneuver, while and `M`
indicates the number of rotor systems.

Every implementation of `AbstractManeuver` must have the properties:

 * `angle::NTuple{N, Function}` where `angle[i](t)` returns the angle of the
        i-th tilting system at time `t` (t is nondimensionalized by the total
        time of the maneuver, from 0 to 1, beginning to end).
 * `RPM::NTuple{M, Function}` where `RPM[i](t)` returns the normalized RPM of
        the i-th rotor system at time `t`. This RPM values are normalized by the
        an arbitrary RPM value (usually RPM in hover or cruise).
"""
abstract type AbstractManeuver{N, M} end

##### FUNCTIONS REQUIRED IN IMPLEMENTATIONS ####################################
"""
    `calc_dV(maneuver::AbstractManeuver, vehicle::Vehicle, t, dt, ttot, Vref)`

Returns the change in velocity `dV=[dVx, dVy, dVz]` (m/s) of `vehicle` performing
`maneuver` at time `t` (s) after a time step `dt` (s).  `Vref` and `tot` are the
reference velocity and the total time at which this maneuver is being performed,
respectively. `dV` is in the global reference system.
"""
function calc_dV(self::AbstractManeuver, vehicle::AbstractVehicle, t::Real,
                                            dt::Real, ttot::Real, Vref::Real)
    error("$(typeof(self)) has no implementation yet!")
end

"""
    `calc_dw(maneuver::AbstractManeuver, vehicle::Vehicle, t, dt, Vref, ttot)`

Returns the change in angular velocity `dW=[dWx, dWy, dWz]` (about global
axes, in radians) of `vehicle` performing `maneuver` at time `t` (s) after a
time step `dt` (s). `ttot` is the total time at which this maneuver is to be
performed.
"""
function calc_dW(self::AbstractManeuver, vehicle::AbstractVehicle, t::Real,
                                                        dt::Real, ttot::Real)
    error("$(typeof(self)) has no implementation yet!")
end

##### COMMON FUNCTIONS  ########################################################
"""
    `get_ntltsys(self::AbstractManeuver)`
Return number of tilting systems.
"""
get_ntltsys(self::AbstractManeuver) = typeof(self).parameters[1]

"""
    `get_nrtrsys(self::AbstractManeuver)`
Return number of rotor systems.
"""
get_nrtrsys(self::AbstractManeuver) = typeof(self).parameters[2]

"""
    `get_angle(maneuver::AbstractManeuver, i::Int, t::Real)`

Returns the angle (in degrees) of the i-th tilting system at the non-dimensional
time t.
"""
function get_angle(self::AbstractManeuver, i::Int, t::Real)
    if i<=0 || i>get_ntltsys(self)
        error("Invalid tilting system #$i (max is $(get_ntltsys(self))).")
    end
    if t<0 || t>1
        @warn("Got non-dimensionalized time $(t).")
    end
    return self.angle[i](t)
end

"""
    `get_angles(maneuver::AbstractManeuver, t::Real)`

Returns the angle (in degrees) of every tilting systems at the non-dimensional
time t.
"""
get_angles(self::AbstractManeuver, t::Real) = Tuple(a(t) for a in self.angle)

"""
    `get_RPM(maneuver::AbstractManeuver, i::Int, t::Real)`

Returns the normalized RPM of the i-th rotor system at the non-dimensional time
t.
"""
function get_RPM(self::AbstractManeuver, i::Int, t::Real)
    if i<=0 || i>get_nrtrsys(self)
        error("Invalid rotor system #$i (max is $(get_nrtrsys(self))).")
    end
    if t<0 || t>1
        @warn("Got non-dimensionalized time $(t).")
    end
    return self.RPM[i](t)
end

"""
    `get_RPMs(maneuver::AbstractManeuver, t::Real)`

Returns the normalized RPM of every rotor systems at the non-dimensional time
t.
"""
get_RPMs(self::AbstractManeuver, t::Real) = Tuple(rpm(t) for rpm in self.RPM)


##### COMMON INTERNAL FUNCTIONS  ###############################################

##### END OF ABSTRACT MANEUVER #################################################










################################################################################
# KINEMATIC MANEUVER TYPE
################################################################################
"""
    `KinematicManeuver(angle, RPM, Vvehicle, anglevehicle)`

A vehicle maneuver where the kinematic are prescribed.

# ARGUMENTS
* `angle::NTuple{N, Function}` where `angle[i](t)` returns the angles
        `[Ax, Ay, Az]` (in degrees)of the i-th tilting system at time `t` (t is
        nondimensionalized by the total time of the maneuver, from 0 to 1,
        beginning to end).
* `RPM::NTuple{M, Function}` where `RPM[i](t)` returns the normalized RPM of
        the i-th rotor system at time `t`. This RPM values are normalized by the
        an arbitrary RPM value (usually RPM in hover or cruise).
* `Vvehicle::Function` where `Vvehicle(t)` returns the normalized vehicle
        velocity `[Vx, Vy, Vz]` at the normalized time `t`. Velocity is
        normalized by a reference velocity (typically, cruise velocity).
* `anglevehicle::Function` where `anglevehicle(t)` returns the angles
        `[Ax, Ay, Az]` (in degrees) of the vehicle relative to the global
        coordinate system at the normalized time `t`.
"""
struct KinematicManeuver{N, M} <: AbstractManeuver{N, M}
    angle::NTuple{N, Function}
    RPM::NTuple{M, Function}
    Vvehicle::Function
    anglevehicle::Function
end

# # Implicit N and M constructor
# KinematicManeuver(a::NTuple{N, Function}, b::NTuple{M, Function},
#                     c::Function, d::Function
#                  ) where {N, M} = KinematicManeuver{N, M}(a, b, c, d)


##### FUNCTIONS  ###############################################################
function calc_dV(self::KinematicManeuver, vehicle::AbstractVehicle, t::Real,
                                            dt::Real, ttot::Real, Vref::Real)
    return Vref * (self.Vvehicle((t+dt)/ttot) - self.Vvehicle(t/ttot))
end

function calc_dW(self::KinematicManeuver, vehicle::AbstractVehicle, t::Real,
                                                         dt::Real, ttot::Real)
    prev_W = (self.anglevehicle(t/ttot) - self.anglevehicle((t-dt)/ttot)) / dt
    cur_W = (self.anglevehicle((t+dt)/ttot) - self.anglevehicle(t/ttot)) / dt
    return pi/180 * (cur_W - prev_W)
end


##### INTERNAL FUNCTIONS  ######################################################

##### END OF KINEMATICMANEUVER  ################################################



"""
    MassProp(m, Ixx, Iyy, Izz, Ixz, Ixy, Iyz)

Mass and moments of inertia in the body frame.
Ixx = int(y^2 + z^2, dm)
Ixz = int(xz, dm)

We can model our quadcopter as two thin uniform rods crossed at the origin with a point mass
(motor) at the end of each. With this in mind, it’s clear that the symmetries result in a diagonal
inertia matrix
"""
struct MassProp{TF}
    m::TF
    Ixx::TF
    Iyy::TF
    Izz::TF
end

"""
    ConstantAtmosphere(Wi, Wb, rho, asound, g)

Constant atmospheric properties.
"""
struct Atmosphere{TF,TV}
    Wi::TV
    Wb::TV
    rho::TF
    asound::TF
    g::TF
end


################################################################################
# KINEMATIC MANEUVER TYPE
################################################################################
mutable struct DynamicManeuver{N, M} <: AbstractManeuver{N, M}
    ttot::Real
    tinit::Real
    angle::NTuple{N, Function}
    RPM::NTuple{M, Function}
    RPMref::Real
    mp::MassProp{Float64}
    atm::Atmosphere{Float64,Array{Float64,1}}
    MovementMatrix::Matrix{Float64}
    k::Real
    b::Real
    kd::Real
    l::Real
    dt::Real
    rotors::Int
end

# # Implicit N and M constructor
# DynamicManeuver(a::NTuple{N, Function}, b::NTuple{M, Function}
#                                     ) where {N, M} = DynamicManeuver{N, M}(a, b)


##### FUNCTIONS  ###############################################################
function calc_dV(self::DynamicManeuver, t::Real)
a = Array{Float64}(undef, 3)
thetaprev = Array{Float64}(undef, 3)
xdotprev = Array{Float64}(undef, 3)
omegaprev = Array{Float64}(undef, 3)
for i in 1:4
rpm[i] = self.RPM[i](t)
end
for i in 1:3
xprev[i] = self.MovementMatix[3+i,selfi.nt]
xdotprev[i] = self.MovementMatrix[6+i,selfi.nt]
thetaprev[i] = self.MovementMatrix[9+i,selfi.nt]
omegaprev[i] = self.MovementMatix[12+i,selfi.nt]
end
a = acceleration(rpm, thetaprev, xdotprev, self.atm, self.k, self.mp)
for i in 1:3
self.MovementMatix[i,selfi.nt+1] = a
end
omegadot = Array{Float64}(undef, 3)
I = [self.mp.Ixx   0.0   0.0;
     0.0     self.mp.Iyy 0.0;
     0.0      0.0  self.mp.Izz]
tau = torque(rpm, self.l, self.k, self.b)
omegadot = inv(I) * (tau - cross(omegaprev, I * omegaprev))
omega = Array{Float64}(undef, 3)
omega = omegaprev + dt * omegadot
for i in 1:3
self.MovementMatix[12+i,selfi.nt+1] = omega
end
thetadot = Array{Float64}(undef, 3)
thetadot = omega2thetadot(thetaprev, omega)
theta = Array{Float64}(undef, 3)
theta = thetaprev + dt * thetadot
for i in 1:3
self.MovementMatix[9+i,selfi.nt+1] = theta
end
xdot = Array{Float64}(undef, 3)
xdot = xdotprev + dt * a
x = Array{Float64}(undef, 3)
x = xprev + dt * xdot
for i in 1:3
self.MovementMatix[3+i,selfi.nt+1] = x
end
dv = xdot - xdotprev
return dv
end


function calc_dw(self::DynamicManeuver, t, nt)
a = Array{Float64}(undef, 3)
thetaprev = Array{Float64}(undef, 3)
xdotprev = Array{Float64}(undef, 3)
omegaprev = Array{Float64}(undef, 3)
for i in 1:4
rpm[i] = self.RPM[i](t)
end
for i in 1:3
xprev[i] = self.MovementMatix[3+i,nt]
xdotprev[i] = self.MovementMatrix[6+i,nt]
thetaprev[i] = self.MovementMatrix[9+i,nt]
omegaprev[i] = self.MovementMatix[12+i,nt]
end
a = acceleration(rpm, thetaprev, xdotprev, self.atm, self.k, self.mp)
for i in 1:3
self.MovementMatix[i,nt+1] = a
end
omegadot = Array{Float64}(undef, 3)
I = [self.mp.Ixx   0.0   0.0;
     0.0     self.mp.Iyy 0.0;
     0.0      0.0  self.mp.Izz]
tau = torque(rpm, self.l, self.k, self.b)
omegadot = inv(I) * (tau - cross(omegaprev, I * omegaprev))
omega = Array{Float64}(undef, 3)
omega = omegaprev + dt * omegadot
for i in 1:3
self.MovementMatix[12+i,nt+1] = omega
end
thetadot = Array{Float64}(undef, 3)
thetadot = omega2thetadot(thetaprev, omega)
theta = Array{Float64}(undef, 3)
theta = thetaprev + dt * thetadot
for i in 1:3
self.MovementMatix[9+i,nt+1] = theta
end
xdot = Array{Float64}(undef, 3)
xdot = xdotprev + dt * a
x = Array{Float64}(undef, 3)
x = xprev + dt * xdot
for i in 1:3
self.MovementMatix[3+i,nt+1] = x
end
return  omega - omegaprev
end


##### INTERNAL FUNCTIONS  ######################################################
function bodytoinertial(theta)

    Ri = Array{eltype(theta)}(undef,3, 3)

    cphi, ctht, cpsi = cos.([theta[1], theta[2], theta[3]])
    sphi, stht, spsi = sin.([theta[1], theta[2], theta[3]])

    Ri[1, 1] = cphi*cpsi - ctht*sphi*spsi
    Ri[1, 2] = -cpsi*sphi - cphi*ctht*spsi
    Ri[1, 3] = stht*spsi

    Ri[2, 1] = ctht*sphi*cpsi + cphi*spsi
    Ri[2, 2] = cphi*ctht*cpsi - sphi*spsi
    Ri[2, 3] = -cpsi*stht

    Ri[3, 1] = sphi*stht
    Ri[3, 2] = cphi*stht
    Ri[3, 3] = ctht

    return Ri

end

"""
 Total thrust, giving the square rotation velocity and the coefficient k
"""
function thrust(rpm, k)
sum=0
for i in 1:self.rotors
    sum  = sum + rpm[i]
end
    T = [0.0; 0.0; k*sum]
    return T
end


"""
Torque, giving the square rotation velocity, the distance of rotors from the center
and the coefficients k and b
"""
function torque(rpm, l, k, b)
    tau = [l*k*(rpm[1]-rpm[3]); l*k*(rpm[2]-rpm[4]); b*(rpm[1]-rpm[2]+rpm[3]-rpm[4])]
    return tau
end

"""
Computation of the linear acceleration
"""
function acceleration(rpm, theta, xdot, atm, k, mp)
    gravity= [0.0; 0.0; -atm.g]
    R = bodytoinertial(theta)
    Tb = thrust(rpm, k)
    Ti = R * Tb
    Fd = -kd * xdot
    a = gravity + 1/mp.m  * Ti + 1/mp.m  * Fd
return a
end

"""
    UniformGravitationalField()

Assumes center of mass and center of gravity are coincident.
"""
#struct UniformGravitationalField <: AbstractInertialModel end

function omega2thetadot(theta,omega)
    ct, cp = cos.([theta[2], theta[1]])
    st, sp = sin.([theta[2], theta[1]])

    C = [1.0  0.0  -st;
        0.0   cp  ct*sp
        0.0  -sp  ct*cp]

    thetadot = inv(C)*omega
    return thetadot
end
##### END OF KINEMATICMANEUVER  ################################################

#
