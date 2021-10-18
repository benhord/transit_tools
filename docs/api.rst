.. _api:

API
===
.. automodule:: transit_tools
    :members:

Transit Tools is centered around the `lightcurve` constructor object with many different analysis tools integrated as methods into the class.

Constructor
~~~~~~~~~~~

It is possible to initialize a `lightcurve` object by providing a `lightkurve.LightCurve` object, providing the time, flux, and optionally flux error arrays, or by providing a target name to instruct the constructor to retrieve or extract a light curve from archival data. `lightcurve` also has the capability to simulate a light curve using `BATMAN` and optionally inject it into a specified existing light curve.
		
.. autoclass:: lightcurve
    :members:

Standalone Analysis Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many of the analysis functions that are wrapped into the `lightcurve` object can also be run as standalone functions through the Transit Tools package.

Light Curve Fetching
--------------------

Transit Search
--------------

.. autofunction:: tls_search

Plotting
--------

.. autofunction:: get_starfield

Light Curve Simulation
----------------------

Signal Vetting/Validation
-------------------------

.. autofunction:: calcfpp_tr
       
Utility Functions
~~~~~~~~~~~~~~~~~

There are some functions included in Transit Tools that do not perform analysis of a light curve but instead are simply used as helper functions.

.. autofunction:: catalog_info

.. autofunction:: tessobs_info

.. autofunction:: coord_to_tic

.. autofunction:: known_pls
