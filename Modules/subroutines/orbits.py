import numpy as np

def length(v):
    return np.sqrt(v[0]**2+v[1]**2+v[2]**2)

# a projected onto b (vectors), output vector
def project(a, b):
    unit_b = unit_vector(b)
    return np.dot(a, unit_b)*unit_b

# all scalar
def find_a(v, r, grav_const, mass):
    return 1/( (2/r) - ((v**2)/(grav_const*mass)) )

# all scalar
def find_e(l, a, grav_const, mass):
    x = (l**2)/(grav_const*mass*a)
    return np.sqrt(1-x)
def unit_vector(vector):
    return vector/np.linalg.norm(vector)

# input scalars, output scalar
def find_phi_0(phi, v_r, a, e, r, l):
    sin = (v_r * a * (1-e**2))/(l*e)
    cos = (a * (1-e**2) - r) / (e*r)
    return phi - np.arctan2(sin, cos)

def find_r(a, e, phi_0, steps):
    if e < 1:
        phis = np.linspace(0,2*np.pi, steps)
    else:
        phi_crit = np.arccos(-1/e)+phi_0
        phis = np.linspace(-phi_crit+0.00001, phi_crit-0.00001, steps)
    r = (a*(1-e**2))/(1+e*np.cos(phis-phi_0))
    r_vector  = np.array([r*np.cos(phis-phi_0), r*np.sin(phis-phi_0),np.zeros_like(r)])
    return r_vector 

# particle velocity relative to gal velocity!  we are in center of mass frame!  (so gal pos is at 0,0,0)
def find_orbit(part_pos, part_vel, gal_mass, grav_const, steps=100):
    rad_vector = -part_pos
    v_r = project(part_vel, rad_vector)
    v_t = part_vel-v_r
    v_r_abs = length(v_r)
    
    v = length(part_vel)
    r = length(rad_vector)
    l_vect = np.cross(rad_vector, part_vel)
    l = length(l_vect)
    print('v_r: ',v_r,' v_t: ',v_t)
    print('v: ',v,' r: ',r,' m: ', gal_mass)
    a = find_a(v, r, grav_const, gal_mass)
    e = find_e(l, a, grav_const, gal_mass)
    impact_parameter = a*(1-e)
    print('a: ', a, ' e: ', e, ' impact parameter: ', a*(1-e))
    
    # now project everything explicitly onto the orbital plane    
    
    e_1 = unit_vector(rad_vector)
    print('pos: ',part_pos,' e_1: ',e_1)
    print(unit_vector(v_r))
    e_3 = unit_vector(l_vect)
    e_2 = np.cross(e_1, e_3)
    print('vectors: ',length(e_1),', ',length(e_2),', ',length(e_3))
    print('e_2: ',e_2,', and unit v_t: ',unit_vector(v_t))
    #e_1 = unit_vector(v_r)
    #e_2 = unit_vector(v_t)
    part_pos_plane = np.array([np.dot(e_1, part_pos),
                               np.dot(e_2, part_pos),
                               0
                              ])
    print('positions: ', part_pos_plane)
    
    phi = np.pi
    
    phi_0 = find_phi_0(phi, v_r_abs, a, e, r, l)
    print('phi0: ', phi_0)
    
    r_vec = find_r(a, e, phi_0, steps)
    print('cross: ', (np.dot(np.cross(e_2,e_1),e_3)))
    # now project everything back to 3D
    e_1_3D = np.array([e_1[0], e_2[0], e_3[0]])
    e_2_3D = np.array([e_1[1], e_2[1], e_3[1]])
    e_3_3D = np.array([e_1[2], e_2[2], e_3[2]])
    r_vec_3D = np.array([np.dot(e_1_3D, r_vec),
                         np.dot(e_2_3D, r_vec),
                         np.dot(e_3_3D, r_vec)
                        ])
    return r_vec_3D, impact_parameter