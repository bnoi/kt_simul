Theory behind the model
=======================

.. warning::
    This document is a mix between french and english. It's currently being
    rewritten to be more intelligible.

.. todo::
    AJOUTER INTRODUCTION TRADUIRE EN ANGLAIS JAI FOUTU EN VRAC LES NOTES
    QUE JAI PRIS

Improving the model
-------------------

* Lenght dependance
    fait varier les sites d'attachement ou la force

* améliorer gradient de aurora en ajoutant un paramètre, actuellement :math:`\frac{d\alpha}{d}`

* :math:`k_{\alpha}` (coef d'attachement) dépend aussi du kt opposé

Anisotropie
-----------

* graphe de :math:`\Delta G` en fonction de la position x
* pente anisotropique(pente du graphe) < 0 qd force de déplacement > à la charge (force de friction de ce qui est transporté par la kinésine, ici le chromatide par exemple)
* energie reçu (chimique -> mécanique): :math:`\Delta G = \Delta G_{ATP->ADP}`

* en cas d'anisotropie nulle = diffusion

* evolution de l'anistropie avec le temps : linéaire et pente négative

Oscillateur
-----------

* :math:`v(t) = \frac{d\theta}{dt}(t)`
* notion d'espace de phase : vitesse en fonction de la position (cercle dont le rayon diminue dans le cas d'un pendule par exemple)
* vecteur d'état : ensemble des positions (1D, 2D, 3D, etc)

* notion de régime amorti et suramorti

KT model
--------

* piston : force % à la vitesse
    :math:`\vec{F} = \mu\vec{v}` avec :math:`\mu` : coef de friction
* ressort : force en fonction de la position (etirement)
    :math:`\vec{F} = k(x - x_{0})`

* kinetochores modélisés par un système piston + ressort

Bilan des forces
^^^^^^^^^^^^^^^^

.. image:: ./_static/images/schema.png
    :width: 700px
    :align: center
    :alt: schematic views

Forces at the right SPB
"""""""""""""""""""""""

* Friction forces (viscous drag):  :math:`F_s^f = -\mu_s \dot{x_s}^R`
* Midzone force generators (applied at the right SPB):

    .. math::
        F_{mid} = F_z\left(1 - (\dot{x}^R_s - \dot{x}_s^L)/V_z\right) =
        F_z\left(1 - 2\dot{x}^R_s / V_z\right)

* Total kinetochore microtubules force generators:

    .. math::

        F_{kMT}^T = \sum_{n = 1}^{N}\sum_{m = 1}^{M_k} & - \rho_{nm}^A\,F_k\left( 1 -
          (\dot{x}^A_{nm} - \dot{x}^R_s)/V_k\right)\\
        & + \lambda_{nm}^A\,F_k\left(1 -
          (\dot{x}^A_{nm} + \dot{x}^R_s)/V_k\right)\\
        & - \rho_{nm}^B\,F_k\left( 1 -
          (\dot{x}^B_{nm} - \dot{x}^R_s)/V_k\right)\\
        & + \lambda_{nm}^A\,F_k\left(1 -
          (\dot{x}^B_{nm} + \dot{x}^R_s)/V_k\right)


Forces at the left SPB
""""""""""""""""""""""

Because of the reference frame definition, :math:`\dot{x_s}^R =
-\dot{x_s}^L\,\forall t`. Here we substituted :math:`x_s^L` with :math:`-x_s^R`

* Friction forces (viscous drag):  :math:`F_f^L = \mu_s \dot{x_s}^R`

* Midzone force generators:
    .. math::
        F_{mid}^L = - F_z\left(1 - 2\dot{x}^R_s / V_z\right)

* Total kinetochore microtubules force generators:
    .. math::
        F_{kMT}^T = \sum_{n = 1}^{N}\sum_{m = 1}^{M_k} & - \lambda_{nm}^A\,F_k\left(1 +
          (\dot{x}^A_{nm} + \dot{x}^R_s)/V_k\right)\\
        & - \lambda_{nm}^B\,F_k\left(1 +
          (\dot{x}^B_{nm} + \dot{x}^R_s)/V_k\right)

Forces at centromere :math:`An`
"""""""""""""""""""""""""""""""

* Drag: :math:`F_c^f = -\mu_c \dot{x_n}^A`
* Cohesin bond (Hook spring) restoring force exerted by
    centromere. We want the centromeres to be able to cross each
    over. In one dimension, this introduces a discontinuity. In the
    previous version, the 'swap' mechanism was solving this directly
    (as :math:`x_A` and :math:`x_B` are exchanged). This is not possible any more, as the 'swap' mechanism is now irrelevant, as there is no prefered
    side for a given centromere.}:

    .. math::

        F_{BA} =
        \kappa_c (x_n^B - x_n^A - d_0) &\mathrm{if}\quad x_n^A \leq x_n^B\\
        \kappa_c (x_n^B - x_n^A + d_0) &\mathrm{if}\quad x_n^A > x_n^B\\

  With :math:`F_{AB} = - F_{BA}`.

* Total visco-elastic bond between the centromere A and the attachment
  sites:

    .. math::
        F_v^T = \sum_{m = 1}^{M_k} -\kappa_k(x_n^A - x_{nm}^A)
        - \mu_k(\dot{x}_n^A - \dot{x}_{nm}^A)

Forces at attachment site :math:`Anm`
"""""""""""""""""""""""""""""""""""""

* Visco-elastic bond between the centromere A and the
  attachment sites:

  .. math::
    F_v =  \kappa_k(x_n^A - x_{nm}^A)
    + \mu_k(\dot{x}_n^A - \dot{x}_{nm}^A)

* Kinetochore microtubules force generators:

    .. math::
        F_{kMT}^A &= F_{kMT}^{RA} + F_{kMT}^{LA}\\
        F_{kMT}^{RA} &= \rho_{nm}^A\,F_k\left(1 - \frac{\dot{x}^A_{nm} -
          \dot{x}^R_s}{V_k}\right)\\
        F_{kMT}^{LA} &=  \lambda_{nm}^A\,F_k\left(-1 - \frac{\dot{x}^A_{nm} -
          \dot{x}^L_s}{V_k}\right)\\

With :math:`F_k = 1` and :math:`V_k = 1` (for now on, we are taking :math:`F_k` as
unit force and :math:`V_k` as unit speed), this gives:

    .. math::
        F_{kMT}^A = \rho_{nm}^A\,\left(\dot{x}^R_s - \dot{x}^A_{nm} + 1\right)%
         - \lambda_{nm}^A\,\left(\dot{x}^R_s + \dot{x}^A_{nm} + 1\right)

Eventually, substituting :math:`\lambda^A_{nm} - \rho^A_{nm}` with :math:`\pi_{nm}^A` and :math:`\lambda^A_{nm} + \rho^A_{nm}` with :math:`|\pi_{nm}^A|`:

    .. math::
        F_{kMT}^A =  \pi_{nm}^A(\dot{x}^R_s + 1) - |\pi_{nm}^A|\dot{x}^A_{nm}

Coefficient d'attachement et détachement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The attachment sites attach or detach stochastically with rates :math:`k_a^{R/L}`
and :math:`k_d`, i.e:

    .. math::
      p_{nm}^A = 1 \xrightarrow{\quad k_d \quad} p_{nm}^A = 0 \xrightarrow{\quad k_a^R \quad} p_{nm}^A = 1\\
      p_{nm}^A= -1 \xrightarrow{\quad k_d \quad} p_{nm}^A = 0 \xrightarrow{\quad k_a^L \quad} p_{nm}^A = -1\\

The detachment rate depends on the position of the attached site with
respect to the chromosome center:

  .. math::

    k_d = k_ad_\alpha / d, \mbox{with } d = \left| x^A_{nm} -
    \left(x^A_{n}+ x^B_{n}\right) / 2 \right|

The attachment rate depends on the state of the other attachment
sites:

    .. math::

        k_a^R = k_a\left( 1/2 + \beta\frac{N_n^{AR} - N_n^{AL}}
        {2\left(N_n^{AR} + N_n^{AL}\right)}\right)

In the discrete time step model, the rates are calculated at each
time step for each attachment site.

Parameters and measures
^^^^^^^^^^^^^^^^^^^^^^^

.. image:: ./_static/images/parameters.png
    :width: 700px
    :align: center
    :alt: parameters

.. image:: ./_static/images/measures.png
    :align: center
    :alt: measures


Equation differientielle du premier ordre
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the viscous nucleoplasm, inertia is negligible. Newton first
principle thus reduces to: $ \sum F = 0 $. This force balance equation
can be written for each elements of the spindle.
To simplify further, the equations for the right and left SPBs can be
combined:

    .. math::
        \mu_s\dot{x}^R_s + F_{z}\left(1 - 2\dot{x}^R_s/V_z\right)
        + \sum_{n,m} - \rho_{nm}^A\,\left(\dot{x}^R_s - \dot{x}^A_{nm} +
          1\right) &= 0 \, \mbox{for the right SPB}\\
        \mu_s\dot{x}^R_s - F_{z}\left(1 - 2\dot{x}^R_s/V_z\right)%
        + \sum_{n,m} \lambda_{nm}^A\,\left(\dot{x}^R_s + \dot{x}^A_{nm} +
          1\right) &= 0 \, \mbox{for the left SPB}\\

The difference of those two expressions gives, with the same substitutions as before:

    .. math::
        \mu_s\dot{x}^R_s + 2F_{z}\left(1 - 2\dot{x}^R_s/V_z\right)
        + \sum_{n,m}- (|\pi_{nm}^A|  + |\pi_{nm}^B|)(\dot{x}^R_s + 1)
        + \pi_{nm}^A \dot{x}_{nm}^A + \pi_{nm}^B \dot{x}_{nm}^B= 0

All the equations are gathered together in the system of equations:

    .. math::
        \mathbf{A}\dot{X} + \mathbf{B}X + C = 0

The vector :math:`X` has :math:`1 + 2N(M_k + 1)` elements and is defined as
follow (Note that the left SPB is omitted in :math:`X`.):

    .. math::
        X = \{x_s^R, \{x_n^A, \{x_{nm}^A\},  x_n^B,
        \{x_{nm}^B \}\}\}\mbox{ with } n \in 1 \cdots N
        \mbox{ and } m \in 1 \cdots M_k


In matrix form, we have:

    .. math::

        X = &
        \begin{pmatrix}
          x_s^R\\
          x_n^A\\
          x_{nm}^A\\
          x_n^B\\
          x_{nm}^B\\
        \end{pmatrix} =
        \begin{pmatrix}
          \text{right SPB}\\
          \text{centromere }A, n\\
          \text{attachment site }A, n,m\\
          \text{centromere }B, n\\
          \text{attachment site }B, n,m\\
        \end{pmatrix}\\
        A = &
        \begin{pmatrix}
          - 2 \mu_s - 4 F_z/V_z - \sum (|\pi_{nm}^A| + |\pi_{nm}^B|)& \hdots & \pi_{nm}^A &
          \hdots &  \pi_{nm}^B\\
          \hdots &  -\mu_c - M_k \mu_k& \mu_k & \hdotsfor{2}\\
          \pi_{nm}^A & \mu_k & - \mu_k - |\pi_{nm}^A| & \hdotsfor{2}\\
          \hdotsfor{3} & -\mu_c - M_k \mu_k & \mu_k\\
          \pi_{nm}^B & \hdotsfor{2} & \mu_k & - \mu_k - |\pi_{nm}^B| \\
        \end{pmatrix}, \\
        B = &
        \begin{pmatrix}
          \,0\, & \hdotsfor{4}\\
          \hdots & - \kappa_c - M_k \kappa_k & \kappa_k &
          \kappa_c & \hdots \\
          \hdots & \kappa_k & -\kappa_k &  \hdotsfor{2}\\
          \hdots & \kappa_c & \hdots &
          -\kappa_c - M_k \kappa_k & \kappa_k \\
          \hdotsfor{3}  & \kappa_k & - \kappa_k\\
        \end{pmatrix}\\
        C = &
        \begin{pmatrix}
          2Fz - \sum_{n,m}(|\pi_{nm}^A| + |\pi_{nm}^B|) \\
          - \delta_n \kappa_c d_0\\
          \pi_{nm}^A\\
          \delta_n \kappa_c d_0\\
          \pi_{nm}^B\\
        \end{pmatrix}
        \mathrm{with}\, \delta_n =
        \begin{cases}
          1  &\mathrm{if}\quad  x_n^A < x_n^B\\
          -1 &\mathrm{if}\quad  x_n^A > x_n^B\\
        \end{cases}

As is actually done in the python implementation,
:math:`A` can be decomposed into a time invariant part :math:`A_0` and a
variable part :math:`A_t` with:\\

    .. math::
        A_0 &=
        \begin{pmatrix}
          - 2 \mu_s - 4 F_z/V_z & \hdotsfor{4}\\
          \hdots &  -\mu_c - M_k \mu_k& \mu_k & \hdotsfor{2}\\
          \hdots & \mu_k & - \mu_k & \hdotsfor{2}\\
          \hdotsfor{3} & -\mu_c - M_k \mu_k & \mu_k\\
          \hdotsfor{3} & \mu_k & - \mu_k\\
        \end{pmatrix}\\
        A_t &=
        \begin{pmatrix}
          - \sum (|\pi_{nm}^A| + |\pi_{nm}^B|)& \hdots & \pi_{nm}^A &
          \hdots &  \pi_{nm}^B\\
          \hdotsfor{5}\\
          \pi_{nm}^A & \hdots & - |\pi_{nm}^A| & \hdotsfor{2}\\
          \hdotsfor{5}\\
          \pi_{nm}^B & \hdotsfor{3} & - |\pi_{nm}^B| \\
        \end{pmatrix}\\

For the sake of clarity, :math:`B` can be decomposed in a kinetochore and a
cohesin part, :math:`B = B_c + B_k`:

    .. math::
        B = \kappa_k
        \begin{pmatrix}
            \,0\, & \hdotsfor{4}\\
            \hdots &  - M_k  & 1 & \hdotsfor{2} \\
            \hdots & 1 & -1 &  \hdotsfor{2}\\
            \hdots &  \hdotsfor{2} & - M_k  & 1 \\
            \hdotsfor{3}  & 1 & - 1\\
        \end{pmatrix}
        + \kappa_c
        \begin{pmatrix}
            \,0\, & \hdotsfor{4}\\
            \hdots & - 1 & \hdots & 1  & \hdots \\
            \hdotsfor{5}\\
            \hdots & 1 & \hdots & -1 & \hdots \\
            \hdotsfor{5}\\
        \end{pmatrix}
