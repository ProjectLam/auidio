# About

This is a godot module meant to detect notes. It creates a AudioEffect which you can add to the audio bus. current implemented algorithms are :

- Autocorrelation

# How to use

This is a godot module. meaning you either need to copy the folder in the modules folder of godot engine or use an external folder and use the ```custom_modules``` property in the ```scons``` command during compilation.

**Example**

Folder structure :

```

┣ godot
┣ custom_modules    ┓
                    ┣ auidio
```
Compilation command :

```scons custom_modules="../custom_modules"```
