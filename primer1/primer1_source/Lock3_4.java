public class Lock3_4 {
    private String query = "";
    private boolean available = false;

    public synchronized String get() {
        while (available == false) {
            try {
                wait();
            } catch (InterruptedException e) { }
        }
        available = false;
        notifyAll();
        return query;
    }

    public synchronized void put(String value) {
        while (available == true) {
            try {
                wait();
            } catch (InterruptedException e) { }
        }
        query = value;
        available = true;
        notifyAll();
    }
}
