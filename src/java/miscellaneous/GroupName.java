package miscellaneous;

import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

public class GroupName {
    public static String getGroupName(){
        DateTimeFormatter format = DateTimeFormatter.ofPattern("yyyyMMDDHHmmss");
        LocalDateTime now = LocalDateTime.now();

        return now.format(format);
    }
}
