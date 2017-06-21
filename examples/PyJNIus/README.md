# PyJNIus Example
[**PyJNIus**](https://github.com/kivy/pyjnius) is a Python package that uses Java reflection and JNI to integrate Java classes into Python. The Jupyter Notebook in this directory demonstrates how to load and use the `SeasonalTrendLoess` classes from Python using **PyJNIus**.

First you'll need to install the `Cython` and `pyjnius` packages. See the **PyJNIus** docs for details. 

To run, start Jupyter Notebook in this directory

```bash
    $ jupyter notebook --port 8989 .
```
A home page should pop up in your browser. Load and run the notebook from there.

On the Mac, if you get an error from `import jnius`, possibly followed by a pop-up telling you to install the legacy Java SE 6 runtime, check the Jupyter Notebook logs. If you see:

```
    No Java runtime present, requesting install
```

then see [this pyjnius issue on how to enable the use of JNI with your JVM](https://github.com/kivy/pyjnius/issues/277).
