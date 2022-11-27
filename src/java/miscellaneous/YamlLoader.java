package miscellaneous;

import org.yaml.snakeyaml.Yaml;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Map;

public class YamlLoader {
    public static Map<String, Object> loadConfigFromYaml(String path) throws IOException {
        InputStream inputStream = new FileInputStream(new File(path));

        Yaml yaml = new Yaml();
        return yaml.load(inputStream);
    }
}
