Developer's Guide
=================

CloudnetPy is hosted by Finnish Meteorological Institute (FMI) and
will be used to process cloud remote sensing data in the
ACTRIS research infrastructure. We are happy to welcome the cloud remote sensing community
to provide improvements in the methods and their implementations, writing
tests and fixing bugs.

How to contribute
-----------------

Instructions can be found from `CloudnetPy's Github page <https://github.com/actris-cloudnet/cloudnetpy/blob/master/CONTRIBUTING.md>`_.

Testing
-------

To run the CloudnetPy test suite, first
clone the whole repository from `GitHub
<https://github.com/actris-cloudnet/cloudnetpy>`_:

.. code-block:: console

	$ git clone https://github.com/actris-cloudnet/cloudnetpy

Testing environment
...................

Now, create a virtual environment and install pytest and CloudnetPy:

.. code-block:: console

    $ cd cloudnetpy
    $ python3 -m venv venv
    $ source venv/bin/activate
    (venv) $ pip3 install pytest .

Unit tests
..........

.. code-block:: console

    (venv) $ pytest

End-to-end test
...............

.. code-block:: console

    (venv) $ python3 tetsts/e2e_test.py

.. note::

   Cloudnetpy performs relatively complicated scientific processing, converting
   noisy measurement data into higher level products. Most of the
   Cloudnetpy's low-level functions are unit tested, but it is
   difficult to write unambiguous tests for the high-level API calls.
   However, the quality of the processed files can be at least roughly
   checked using CloudnetPy's quality control functions.


Coding guidelines
-----------------

- Use `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ standard.

- Write `Google-style docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_.

- Check your code using, e.g., `Pylint <https://www.pylint.org/>`_.
