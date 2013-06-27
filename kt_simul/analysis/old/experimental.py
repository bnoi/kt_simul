"""
Initially in simul_spindle.py
"""

# def _make_movie(self, t, imsize):
#     """somehow deprecated"""
#     try:
#         os.mkdir("movie")
#     except OSError:
#         pass

#     imname = "movie/spindle%03i.tif" % t
#     if os.path.isfile(imname):
#         os.remove(imname)

#     imarray = np.zeros(imsize, np.uint8)

#     xspbR = scale(self.KD.spbR.pos, imsize[1])
#     yspbR = scale(0, imsize[0])
#     xspbL = scale(self.KD.spbL.pos, imsize[1])
#     yspbL = scale(0, imsize[0])
#     imarray[yspbR, xspbR, 0] += 128
#     imarray[yspbL, xspbL, 0] += 128

#     chromosomes = self.KD.chromosomes
#     for n, ch in enumerate(chromosomes):
#         xchR = scale(ch.cen_A.pos, imsize[1])
#         ychR = scale(0.1 * n - 0.1, imsize[0])
#         xchL = scale(ch.cen_B.pos, imsize[1])
#         ychL = scale(0.1 * n - 0.1, imsize[0])
#         imarray[ychR, xchR, 1] += 42
#         imarray[ychL, xchL, 1] += 42

#     if 80 < t < 100:
#         self.testim =  imarray
#         print xchR, ychR, t
#     im = image_fromarray(imarray, 'RGB')
#     im.save(imname)
#     del im

# def get_2Dtraj_list(self, ang_noise=5e-3, pos_noise=1e-2,
#                     perp_distance=0.3,
#                     cen=None, cen_sigma=0.6):

#     """
#     Returns a list of 2D trajectories, adding a random mouvement
#     to the whole spindle.

#     Keyword arguments:
#     ------------------
#     ang_noise : float
#         angular variation in radians per seconds
#     pos_noise : float
#         center of mass displacement in um/s
#     perp_distance: float
#         perpendicular distance between the kinetochore pairs
#     cen: int or None
#         chromosome index for which a centromere trajectory is
#         simulated. Be aware that this adds a random coiled coil
#         movement around the kinetochore position.
#         If cen is None returns all 6 trajectories plus the SPBs
#     cen_sigma: float
#         sets the amplitude of the simulated random coiled coil
#         movement of the cen marker with respect to the chromosome

#     Returns:
#     --------
#     trajectories: ndarray
#         Trajectories is an array of shape (n, 2, num_steps),
#         where n is the number of trajectories
#     """

#     dt = self.KD.params['dt']
#     N =  self.KD.params['N']
#     ang_noise *= dt # rad.s^-1
#     pos_noise *=  dt # um.s^-1

#     n_traj = 2 + 2 * N if cen is None else 4
#     trajs_shape = (n_traj, 2, self.num_steps)
#     trajectories = np.zeros(trajs_shape)

#     trajectories[0, 0, :] = self.KD.spbR.traj
#     trajectories[1, 0, :] = self.KD.spbL.traj

#     if cen is not None:
#         ch = self.get_ch(cen)
#         trajectories[2, 0, :] = ch.cen_A.traj
#         trajectories[3, 0, :] = ch.cen_B.traj
#         trajectories[2:, ...] += normal(0, scale=cen_sigma,
#                                        size=(2, 2, self.num_steps))
#     else:
#         for n, ch in enumerate(self.KD.chromosomes):
#             trajectories[2 + 2 * n, 0, :] = ch.cen_A.traj
#             trajectories[2 + 2 * n, 1, :] += (1 - n) * perp_distance
#             trajectories[2 + 2 * n + 1, 0, :] = ch.cen_B.traj
#             trajectories[2 + 2 * n + 1, 1, :] += (1 - n) * perp_distance
#     # TODO: block needs to be vectorized GG march 2012
#     xcs = []
#     ycs = []
#     rots = []
#     xd, yd, thetad = (0.,)*3
#     for n in range(self.num_steps):
#         xd += normal(0, scale = pos_noise)
#         yd += normal(0, scale = pos_noise)
#         thetad += normal(0, scale = ang_noise)
#         xcs.append(xd)
#         ycs.append(yd)
#         rots.append([[np.cos(thetad), np.sin(thetad)],
#                      [np.cos(thetad), - np.sin(thetad)]])
#     # TODO: block needs to be vectorized GG march 2012
#     for traj in trajectories:
#         traj += np.vstack((np.array(ycs), np.array(xcs)))
#         n = 0
#         for pos, rot in zip(traj, rots):
#             new_pos = np.dot(pos, rot)
#             traj[n] = new_pos
#             n += 1
#     # TODO: implement a general vectorized Brownian motion
#     return trajectories

# def get_3Dtraj_list(self, ang_noise=5e-2, pos_noise=3e-2,
#                     radial_distance=0.3, cen=None, cen_sigma=0.6):
#     """
#     Returns a list of 3D trajectories, adding a random mouvement [1]_
#     to the whole spindle.

#     Keyword arguments:
#     ------------------
#     ang_noise : float
#         angular variation in radians per seconds, sets the standard
#         deviation of the spindle axis angle in um/s.

#     pos_noise : float
#         std. dev of the spindle center displacement in um/s
#     radial_distance: float
#         radial distance of the centromeres to the spindle axis.
#         Here, chromosomes are distributed evenly around the spindle
#         axis.
#     cen: int or None
#         chromosome index for which a centromere trajectory is
#         simulated. Be aware that this adds a random coiled coil
#         movement around the kinetochore position.
#         If cen is None returns all 6 trajectories plus the SPBs
#     cen_sigma: float or None
#         if cen is not None, sets the

#     Returns:
#     --------

#     trajectories: a list of ndarrays


#     .. [1] By adding a normaly distributed noise with a Gaussian
#            distribution see numpy.random.normal for further details
#     """

#     dt = self.KD.params['dt']
#     N = self.KD.params['N']
#     trajectories = []
#     ang_noise *= dt # rad.s^-1
#     pos_noise *=  dt # um.s^-1

#     n_traj = 2 + 2 * N if cen is None else 4
#     trajs_shape = (n_traj, 3, self.num_steps)
#     trajectories = np.zeros(trajs_shape)
#     trajectories[0, 0, :] = self.KD.spbR.traj
#     trajectories[1, 0, :] = self.KD.spbL.traj

#     if cen is not None:
#         ch = self.get_ch(cen)
#         trajectories[2, 0, :] = ch.cen_A.traj
#         trajectories[3, 0, :] = ch.cen_B.traj
#         trajectories[2:, ...] += normal(0, scale=cen_sigma,
#                                        size=(2, 2, self.num_steps))
#     else:
#         for n, ch in enumerate(self.KD.chromsomes):
#             trajectories[2 + 2 * n, 0, :] = ch.cen_A.traj
#             phi = (n - 1) * 2 * np.pi / N
#             radial_y = radial_distance * np.cos(phi)
#             trajectories[2 + 2 * n, 1, :] += radial_y
#             trajectories[2 + 2 * n + 1, 1, :] += radial_y
#             radial_z = radial_distance * np.sin(phi)
#             trajectories[2 + 2 * n + 1, 2, :] += radial_z

#     xcs, ycs, zcs = np.zeros((3, self.num_steps))
#     #Fix the x axis, rotate around the two others
#     xy_rots = np.zeros((3, 3, self.num_steps))
#     zx_rots = np.zeros((3, 3, self.num_steps))

#     thetad, phid = 0., 0.
#     for n in range():
#         thetad += normal(0, scale = ang_noise)
#         phid += normal(0, scale = ang_noise)
#         xcs[n] += normal(0, scale = pos_noise)
#         ycs[n] += normal(0, scale = pos_noise)
#         zcs[n] += normal(0, scale = pos_noise)
#         xy_rots[:, :, n] = [[np.cos(thetad), - np.sin(thetad), 0],
#                             [np.sin(thetad), np.cos(thetad), 0],
#                             [0, 0, 1]]
#         zx_rots[:, :, n] = [[np.cos(phid), 0, - np.sin(phid)],
#                             [0, 1, 0],
#                             [np.sin(phid), 0, np.cos(phid)]]
#     for traj in trajectories:
#         traj += np.vstack((xcs, ycs, zcs))
#         for n, pos in enumerate(traj.T):
#             tmp_pos = np.dot(pos, xy_rots[:, :, n])
#             new_pos = np.dot(tmp_pos, zx_rots[:, :, n])
#             traj[:, n] = new_pos
#     return trajectories

# def get_ch(self, n = 0):
#     return self.KD.chromosomes[n]

# def scale(x, size, pix_size = 0.0645):
#     """
#     Scale the position x on a line of size "size" from microns to pixels
#     origin is put in the center of the line
#     """
#     return int(x / pix_size + size / 2)