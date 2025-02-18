run("Bio-Formats Macro Extensions");

// Selecciona un directorio de imágenes para abrir
input = getDirectory("Pon un folder donde están las imágenes");

// Selecciona un directorio para las nuevas imágenes
output = getDirectory("Folder de nuevas imágenes");

// Obtiene la lista de archivos en el directorio de entrada
list = getFileList(input);

// Activa el modo por lotes
setBatchMode(true);

// Procesa cada archivo en la lista
for (i = 0; i < list.length; i++) {
    Ext.openImagePlus(input + list[i]); // Abre el archivo actual

    // Aplica las siguientes operaciones de procesamiento
run("8-bit");
run("Grays");
run("Enhance Contrast...", "saturated=0.35 process_all");
run("Unsharp Mask...", "radius=1.2 mask=0.5 stack");
run("Despeckle", "stack");
run("Subtract Background...", "rolling=40 stack");
run("Z Project...", "projection=[Max Intensity]");
// Guarda la imagen procesada en el directorio de salida
    saveAs("Tiff", output + list[i]);
}

// Desactiva el modo por lotes
setBatchMode(false);