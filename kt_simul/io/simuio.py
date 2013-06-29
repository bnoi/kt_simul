"""
Main module for save and read simulation.
"""

import time
import logging
import StringIO
import zipfile

import numpy as np

from xml.etree.ElementTree import Element, SubElement, tostring

from kt_simul.core.spindle_dynamics import KinetoDynamics
from kt_simul.core.simul_spindle import Metaphase
from kt_simul.io.xml_handler import ParamTree, indent, ResultTree

logger = logging.getLogger(__name__)


class SimuIO():
    """
    Parameters
    ----------
    meta_instance : :class:`~kt_simul.core.simul_spindle.Metaphase`
        Should have been already perform a simulation

    """

    def __init__(self, meta_instance=None):

        if meta_instance:
            self.meta = meta_instance
            self.KD = self.meta.KD
            self.timelapse = self.meta.timelapse
            self.num_steps = self.meta.num_steps
            self.paramtree = self.meta.paramtree
            self.measuretree = self.meta.measuretree
            self.observations = self.meta.observations

    def save(self, simufname="results.zip", verbose=False):
        """
        Saves the results of the simulation in two files
        with the parameters, measures and observations in one file
        and the trajectories in the other.

        Parameters
        ----------
        simufname : string
            A ZIP file which contains results.xml and data.npy:
                - results.xml : name of the xml file where parameters and
                                observations will be written

                - data.npy : data are stored in numpy's binary format

        """

        chromosomes = self.KD.chromosomes

        # Store matrices insisde wavelist
        wavelist = []

        xmlout = ""
        xmlout += '<?xml version="1.0"?>\n'
        today = time.asctime()
        experiment = Element("experiment", date=today)
        experiment.append(self.paramtree.root)
        experiment.append(self.measuretree.root)

        # right SPB
        spbR = SubElement(experiment, "trajectory", name="rightspb",
                          column='0', units='mu m')
        SubElement(spbR, "description").text = "right spb trajectory"
        spbRtraj = np.array(self.KD.spbR.traj)
        wavelist.append(spbRtraj)

        # left SPB
        spbL = SubElement(experiment, "trajectory", name="leftspb",
                          column='1', units='mu m')
        SubElement(spbL, "description").text = "left spb trajectory"
        spbLtraj = np.array(self.KD.spbL.traj)
        wavelist.append(spbLtraj)

        col_num = 2
        # chromosomes
        for n, ch in enumerate(chromosomes):
            rch = SubElement(experiment, "trajectory",
                            name="centromereA",
                            index=str(n),
                            column=str(col_num),
                            units='mu m')
            text = "chromosome %i centromere A trajectory" % n
            SubElement(rch, "description").text = text
            wavelist.append(ch.cen_A.traj)
            col_num += 1

            SubElement(experiment, "numbercorrect", name="centromereA",
                       index=str(n), column=str(col_num))
            wavelist.append(ch.correct_history[:, 0])
            col_num += 1
            SubElement(experiment, "numbererroneous",
                       name="centromereA", index=str(n),
                       column=str(col_num))
            wavelist.append(ch.erroneous_history[:, 0])
            col_num += 1

            lch = SubElement(experiment, "trajectory", index=str(n),
                             column=str(col_num), units='mu m')
            text = "chromosome %s left kinetochore trajectory" % n
            SubElement(lch, "description").text = text
            wavelist.append(np.array(ch.cen_B.traj))
            col_num += 1
            SubElement(experiment, "numbercorrect", name="centromereB",
                       index=str(n), column=str(col_num))
            wavelist.append(ch.correct_history[:, 1])
            col_num += 1

            SubElement(experiment, "numbererroneous", name="centromereB",
                       index=str(n), column=str(col_num))
            wavelist.append(ch.erroneous_history[:, 1])
            col_num += 1

            #Plug Sites
            for m, plugsite in enumerate(ch.cen_A.plugsites):
                SubElement(experiment, "trajectory", name="plugsite",
                           index=str((n, m)), cen_tag='A',
                           column=str(col_num), units='mu m')
                wavelist.append(plugsite.traj)
                col_num += 1

                SubElement(experiment, "state", name="plugsite",
                           index=str((n, m)), cen_tag='A',
                           column=str(col_num), units='')
                wavelist.append(plugsite.state_hist)
                col_num += 1

            for m, plugsite in enumerate(ch.cen_B.plugsites):

                SubElement(experiment, "trajectory", name="plugsite",
                           index=str((n, m)), cen_tag='B',
                           column=str(col_num), units='mu m')
                wavelist.append(plugsite.traj)
                col_num += 1

                SubElement(experiment, "state", name="plugsite",
                           index=str((n, m)), cen_tag='B',
                           column=str(col_num), units='')
                wavelist.append(plugsite.state_hist)
                col_num += 1

        # Observations
        obs_elem = SubElement(experiment, "observations")
        if not hasattr(self, 'observations'):
            self.evaluate()
        for key, val in self.observations.items():
            SubElement(obs_elem, key).text = str(val)

        # Now we write down the whole experiment XML element
        indent(experiment)
        xmlout += tostring(experiment)

        # And the data in the file datafname
        datafile = StringIO.StringIO()
        data = np.vstack(wavelist).T
        np.save(datafile, data)

        # Pack xmlfile and datafile in a zip archive
        zipf = zipfile.ZipFile(simufname, "w",
            compression=zipfile.ZIP_DEFLATED)
        zipf.writestr("results.xml", xmlout)
        zipf.writestr("data.npy", datafile.getvalue())
        zipf.close()

        datafile.close()

        if verbose:
            logger.info("Simulation saved to file %s " % simufname)

    def read(self, simufname="results.zip", verbose=False):
        """
        Creates a simul_spindle.Metaphase from a results.zip file.

        Parameters
        ----------
        simufname : string, optional
            The .zip file where results from the existing simulation are (zip which contains xml and npy)

        verbose : bool
            Set Metaphase verbose

        Returns
        -------
        :class:`~kt_simul.core.simul_spindle.Metaphase`

        """

        # Unzip results.xml and data.npy and put it in
        # tempfile
        try:
            zipf = zipfile.ZipFile(simufname, "r")
        except:
            logger.info("%s does not appear to be a Zipfile" % simufname)
            return False

        # Test if simu results file are in the archive
        if "results.xml" not in zipf.namelist() or "data.npy" not in zipf.namelist():
            logger.info("%s does not contain results.xml and data.npy" % simufname)
            return False

        xmltemp = StringIO.StringIO(zipf.read("results.xml"))
        datatemp = StringIO.StringIO(zipf.read("data.npy"))

        restree = ResultTree(xmltemp, datatemp)

        param_root = restree.root.find('parameters')
        paramtree = ParamTree(root=param_root)
        params = paramtree.relative_dic
        measure_root = restree.root.find('measures')
        measuretree = ParamTree(root=measure_root,
                                adimentionalized=False)

        metaphase = Metaphase(paramtree=paramtree,
                              measuretree=measuretree,
                              verbose=verbose)

        traj_matrix = restree.get_all_trajs()
        correct_matrix = restree.get_all_correct()
        erroneous_matrix = restree.get_all_erroneous()
        state_hist_matrix = restree.get_all_plug_state()
        KD = KinetoDynamics(params)
        KD.spbR.traj = traj_matrix[:, 0]
        KD.spbL.traj = traj_matrix[:, 1]
        col_num = 2
        state_num = 0
        for n, ch in enumerate(KD.chromosomes):
            ch.cen_A.traj = traj_matrix[:, col_num]
            ch.cen_A.pos = ch.cen_A.pos
            col_num += 1
            ch.cen_B.traj = traj_matrix[:, col_num]
            ch.cen_B.pos = ch.cen_B.pos
            col_num += 1
            ch.erroneous_history = (erroneous_matrix[:, n * 2: n * 2 + 2])
            ch.correct_history = (correct_matrix[:, n * 2: n * 2 + 2])
            for plugsite in ch.cen_A.plugsites:
                plugsite.traj = traj_matrix[:, col_num]
                col_num += 1
                plugsite.state_hist = state_hist_matrix[:, state_num]
                plugsite.plug_state = plugsite.state_hist[-1]
                state_num += 1
            for plugsite in ch.cen_B.plugsites:
                plugsite.traj = traj_matrix[:, col_num]
                col_num += 1
                plugsite.state_hist = state_hist_matrix[:, state_num]
                plugsite.plug_state = plugsite.state_hist[-1]
                state_num += 1

        metaphase.KD = KD

        metaphase.KD.simulation_done = True

        return metaphase
