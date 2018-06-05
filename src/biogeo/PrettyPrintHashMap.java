package biogeo;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Arrays;

public class PrettyPrintHashMap<K, V> {
	
	private HashMap<K, V> map;
	
	public PrettyPrintHashMap(HashMap<K, V> map) {
		this.map = map;
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		Iterator<Entry<K, V>> iter = map.entrySet().iterator();
		
		while (iter.hasNext()) {
			Entry<K, V> entry = iter.next();
			String[] key = (String[]) entry.getKey();
			sb.append(Arrays.toString(key));
			sb.append(" = ");
			sb.append(entry.getValue());
			
			if (iter.hasNext()) {
				sb.append(',').append(' ');
			}
		}
		
		return sb.toString();
	}
}
