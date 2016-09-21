public class MatrixMapper extends
Mapper<LongWritable, Text, Text, Text>
{
  @Override
  protected void map
  (LongWritable key, Text value, Context context)
  throws IOException, InterruptedException {
  String[] csv = value.toString().split(",");
  String matrix = csv[0].trim();
  int row = Integer.parseInt(csv[1].trim());
  int col = Integer.parseInt(csv[2].trim());
  if(matrix.contains("a")) {
    for (int i=0; i < lMax; i++) {
      String akey = Integer.toString(row) + "," +
        Integer.toString(i);
      context.write(new Text(akey), value);
    }
  }
  if(matrix.contains("b")) {
    for (int i=0; i < iMax; i++) {
      String akey = Integer.toString(i) + "," +
        Integer.toString(col);
      context.write(new Text(akey), value);
    }
  }
  }
}

public class MatrixReducer extends
Reducer<Text, Text, Text, IntWritable> {
  @Override
  protected void reduce
  (Text key, Iterable<Text> values, Context context)
  throws IOException, InterruptedException {
  int[] a = new int[5];
  int[] b = new int[5];

  for (Text value : values) {
    System.out.println(value);
    String cell[] = value.toString().split(",");
    if (cell[0].contains("a")) {
      int col = Integer.parseInt(cell[2].trim());
      a[col] = Integer.parseInt(cell[3].trim());
    }
    else if (cell[0].contains("b")) {
      int row = Integer.parseInt(cell[1].trim());
      b[row] = Integer.parseInt(cell[3].trim());
    }
  }
  int total = 0;
  for (int i = 0; i < 5; i++) {
    int val = a[i] * b[i];
    total += val;
  }
  context.write(key, new IntWritable(total));
  }
}
