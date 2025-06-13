import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

/**
 * Checks if a file is gzipped by reading its first two bytes.
 *
 * @param file The file to check (can be File or Path).
 * @return true if the file is gzipped, false otherwise.
 */
boolean isGzipped(Object file) {
    // Convert to Path for unified handling
    Path path
    if (file instanceof File) {
        path = file.toPath()
    } else if (file instanceof Path) {
        path = file
    } else {
        path = Paths.get(file.toString())
    }
    
    // Check if file exists and has sufficient length
    if (!Files.exists(path) || Files.size(path) < 2) {
        return false
    }
    
    try {
        byte[] header = new byte[2]
        Files.newInputStream(path).withCloseable { stream ->
            stream.read(header, 0, 2)
        }
        // GZIP magic number is 0x1f8b
        return header[0] == 0x1f && header[1] == (byte) 0x8b
    } catch (IOException e) {
        // Handle potential I/O errors
        e.printStackTrace()
        return false
    }
}