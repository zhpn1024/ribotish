=================================
INSTALL Guide For Ribo-TISH
=================================

Please check the following instructions to complete your installation.

Prerequisites
=============

Python version must be equal to *2.7*.

:Pysam_: >= 0.8.3
:Scipy_: >= 0.15.1
:Matplotlib_: >= 1.4.3

.. _Pysam: https://pypi.python.org/pypi/pysam
.. _Scipy: http://www.scipy.org/Download
.. _Matplotlib: http://matplotlib.org/users/installing.html

Ribo-TISH cannot run on Windows because Pysam_ do not support Windows installation.

Easy installation through PyPI
==============================

The easiest way to install Ribo-TISH is through PyPI system. Get pip_ if it's not available in your system. 

Then under command line::

  $ pip install ribotish

PyPI will install all dependents automatically if they are absent.  

To upgrade Ribo-TISH, type::

  $ pip install -U ribotish

If you do not want pip to fix dependencies, type::

  $ pip install --no-deps ribotish

To install/upgrade Ribo-TISH under HOME directory without fixing dependencies::

  $ pip install --no-deps -U --user ribotish

.. _pip: http://www.pip-installer.org/en/latest/installing.html

Install from source
===================

To install a source distribution, download source from PyPI_ or GitHub_. The latest updates usually appear on GitHub first.::

 $ git clone https://github.com/zhpn1024/ribotish

Go to the directory and simply run the install script::

 $ python setup.py install

By default, the script will install python library and executable codes globally, which means you need to be root or administrator of the machine so as to complete the installation. 

If I want to install everything under my own HOME directory, use this command::

 $ python setup.py install --prefix /home/pzhang/

or::

 $ python setup.py install --user

.. _PyPI: https://pypi.python.org/pypi/ribotish
.. _GitHub: https://github.com/zhpn1024/ribotish

Configure enviroment variables
==============================

After running the setup script, you might need to add the install location to your ``PYTHONPATH`` and ``PATH`` environment variables. The process for doing this varies on each platform, but the general concept is the same across platforms.

PYTHONPATH
~~~~~~~~~~

To set up your ``PYTHONPATH`` environment variable, you may need to add the value ``PREFIX/lib/pythonX.Y/site-packages`` to your existing ``PYTHONPATH``. In this value, X.Y stands for the major–minor version of Python you are using (such as 2.7). ``PREFIX`` is the install prefix where you installed ribotish. If you did not specify a prefix on the command line, it will be installed using Python's sys.prefix value.

On Linux, using bash, you can include the new value in ``PYTHONPATH`` by
adding this line to ``~/.bashrc`` or ``~/.bash_profile``. For example::

 $ export PYTHONPATH=~/lib/python2.7/site-packages:$PYTHONPATH


PATH
~~~~

Just like your ``PYTHONPATH``, you may also need to add a new value to your PATH environment variable so that you can use the ribotish command line directly. This time you need to add ``PREFIX/bin`` to your PATH environment variable. The process is the same as for the ``PYTHONPATH`` variable::

 $ export PATH=~/bin:$PATH

