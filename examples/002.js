
let data = {
  labels: ["Mon", "Tue", "Wed", "Thur", "Fri"],
  series: [[1, 3, 5, 3, 7]],
};

let options = {
  width:900,
  height:300,
};

new Chartist.Bar('.ct-chart', data, options);
