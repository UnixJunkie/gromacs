Bugs fixed
^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

The deform option was unsuited for flow simulations
"""""""""""""""""""""""""""""""""""""""""""""""""""

The deform option deformed the coordinates, in addition to the box, and did
not correct the velocities of particles when they were shifted by a periodic
box vector. This has now been corrected which makes the deform option also
useful for shear flows. Applications were the system was strechted until
some interactions broke were probably not affected measurably by
these issues. Note that a velocity profile should be generated when using
deform with the current or later versions. An mdp option has been added
to let ``grompp`` do this.

:issue:`4607`
