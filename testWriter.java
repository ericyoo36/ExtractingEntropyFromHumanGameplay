import java.io.BufferedWriter;
import java.io.FileWriter;
import java.security.SecureRandom;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

public class testWriter {

	public static void main(String[] args) {
		SecureRandom rand = new SecureRandom();
		try {
		BufferedWriter writer = new BufferedWriter(new FileWriter("test.txt"));
		int x = 0;
		int y = 0;
		for (int i = 0; i < 40000; i++) {
			x = rand.nextInt(1000);
			y = rand.nextInt(600);
			//writer.write("X: " + ThreadLocalRandom.current().nextInt(0, 1090) + "\tY: " + ThreadLocalRandom.current().nextInt(0, 580) + "\tTicks: 16\n");
			writer.write("X: " + x + "\tY: " + y + "\tTicks: 16\n");
		}
		
		writer.close();
		
		} catch (Exception e) {
			System.out.println("writeTest Failed");
			return;
		}
	}
}
