package main;

import java.io.File;
import java.util.*;

public class SuperSetRemoval {

	public static void main(String[] args) throws Exception{
		// TODO Auto-generated method stub
		long startTime = System.nanoTime();
		
		ArrayList<Long> temp = new ArrayList<>();
		File f = new File("to_send_to_java.txt");
		
		Scanner fileReader = new Scanner(f);
		int n = 0;
		while(fileReader.hasNextLine()) {
			n++;
			temp.add(Long.parseLong(fileReader.nextLine(), 2));
		}
		Long arr[] = temp.toArray(new Long[n]);
		fileReader.close();
		boolean removed[] = new boolean[n];
		
		int reds = 0;
		int supers = 0;
		int MCSs = 0;		
		for (int i=0; i<n; i++) {
			if((i>>16)<<16==i) {
				System.err.printf("\n%.2f",((float)i/n)*((float)i/n)*100 );
				System.err.println("% completed \nseconds passed: " + (System.nanoTime() - startTime)/1000000000);
				System.err.println(reds);
				System.err.println(supers);
			}
			if (removed[i])
				continue;
			for(int j=0; j<i; j++) {
				if (removed[j])
					continue;
				if (arr[i]==arr[j]) {
					reds++;
					removed[j] = true;
					continue;
				}
				if ((arr[i]&arr[j])==arr[i]) {
					supers++;
					removed[j] = true;
					continue;
				}
				if ((arr[i]&arr[j])==arr[j]) {
					supers++;
					removed[i] = true;
					break;
				}
				
			}
			if (!removed[i])
				MCSs++;
		}
		
		long endTime = System.nanoTime();
		System.err.println("number of repeated MCSs: "+reds);
		System.err.println("number of non-MCSs(Super set of a MCS a.k.a cut-sets)" + supers);
		System.out.println("duration of the post-processing: " + (endTime - startTime)/1000000 + " miliseconds");
		System.out.println("Total MCSs: " + MCSs);
	}

}
