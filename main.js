class Dataset{

    constructor(){
        // data to read
        this.wavels = [];
        this.spectra = [];
        this.headers = [];

        // derived from read data
        this.wavel_max = 0;
        this.wavel_min = 0;
    }

    _getWavelMinMax(){
        // we'll need these for interpolation later
        this.wavel_max = Math.max(...this.wavels);
        this.wavel_min = Math.min(...this.wavels);
    }

    getDyeVector(dye_index){
        // return an object that tells you which dye vector is returned and then also the data
        let dye_data = this.spectra.map(row => row[dye_index]);
        let dye_name = this.headers[dye_index];

        // returns an object as dye_name and dye_data
        return {dye_name, dye_data}
    }

    loadFromContents(csv_contents){

        // split the str into lines
        let lines = csv_contents.trim().split('\n');

        // get the numer of lines i need to iterate over to read the dat
        const num_lines = lines.length

        for (let nth_line = 0; nth_line < num_lines; nth_line++){
                
            // get the nth line 
            let line = lines[nth_line].split(',');

            // grab the data
            if (nth_line > 0){
                // add the wavelength value to the wavel vector
                this.wavels.push(+line[0]);
                // add the data vector for each of the dyes to the nth row of the data array
                // essentially represeting a 2D array as a list of lists
                this.spectra.push(line.slice(1).map(Number));
            }
            // grab the dye names
            else if (nth_line == 0){
                // this will add the array of sample names to the headers attribute
                this.headers = line.slice(1);
            }
        }

        // get the wavelength range ends 
        this._getWavelMinMax()

        return this
    }

    loadFromVariables(wavels, spectra, headers){

        this.wavels = wavels
        this.spectra = spectra
        this.headers = headers
        this._getWavelMinMax()
        
        return this
    }
}

class Interpolator{

    constructor(dataset, query_points){
        // theres probably a better way to do this than passing in a dataset obj
        // and returning a interpolated dataset

        // get the two x vectors we'll need
        this.input_wavels = dataset.wavels;
        this.query_points = query_points;
        this.headers = dataset.headers

        // this will be an m by n array 
        // where m is the num_wavels and n is num_dyes
        // each dye spectra will have to be interpolated
        this.spectra = dataset.spectra;
        
        // also would be good to have the number of dyes
        this.num_dyes = dataset.spectra[0].length;
        
        //-----------------------------------------------------------------
        // Do the interpolation here
        this.interpolated_spectra = []
        this.interpolateDataset()

        //-----------------------------------------------------------------

    }

    _linearSearch(x_q){
        // takes some query point in question (x_q) and searches the input wavels 
        // for the indices x_i and x_i+1 where x_i <= x_q <= x_i+1
        for (let ind = 0; ind < this.input_wavels.length - 1; ind++){
            
            // here here we can return the values x_i and x_i+1
            if ((this.input_wavels[ind] <= x_q) & (x_q <= this.input_wavels[ind + 1])){
                
                // going to prepare an output vector here with the things we need to interpolate
                // basically; [xi, Yi, xi+1, Yi+1, xq]
                // in this case Y is a vector corresponding to the y value for the nth dye
                let output_vector = [this.input_wavels[ind],this.spectra[ind],this.input_wavels[ind+1],this.spectra[ind+1],x_q];

                return output_vector
            }
        } 
    }

    _interpolatePoints(output_vector){

        // just going to unpack the output vector real quick
        let xi = output_vector[0];
        let Yi = output_vector[1];
        let xi_1 = output_vector[2];
        let Yi_1 = output_vector[3];
        let xq = output_vector[4];
        
        // make a vector to hold the n interpolated points at wavelength xq and m and yq
        let Yq = [];
        let m = 0;
        let yq = 0;

        // so the Y vectors are coming in with n elements where n corresponds to the nth dye 
        for (let nth_dye = 0; nth_dye < this.num_dyes; nth_dye++){

            // so for linear interpolation we'll use:
            // yq = m(xq-xi)+yi 
            // where m = (yi+1 - yi) / (xi+1 - xi)
            m = (Yi_1[nth_dye] - Yi[nth_dye]) / (xi_1 - xi);
            
            // for the nth dye the interpolated emission is xq
            yq = m * (xq - xi) + Yi[nth_dye];
            
            // add the nth_result
            Yq.push(yq);
        }
        return Yq
    }

    interpolateDataset(){

        // need to reset the results so its not getting doubled and stuff
        this.interpolated_spectra = []
        
        // init
        let xq = 0;
        let search_results = [];
        let Yq = [];

        // interpolate each query_point 
        for (let ind = 0; ind < this.query_points.length; ind++){

            // get xq 
            xq = this.query_points[ind];
            // find xi and xi+1
            search_results = this._linearSearch(xq);
            // do linear interpolation between these points
            Yq = this._interpolatePoints(search_results);
            // add the newly interpolated point to the interpolated spectra set
            this.interpolated_spectra.push(Yq);
        }
        let interpolated_dataset = new Dataset()
        
        return interpolated_dataset.loadFromVariables(this.query_points, this.interpolated_spectra, this.headers)

    }
} 

function readFile(file) {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (e) => resolve(e.target.result);
        reader.onerror = () => reject(new Error("Error reading file."));
        reader.readAsText(file);
    });
}

class Matrix {

    constructor(data){
        // only supporting 2D matricies
        this.data = data
        this.num_rows = data.length;
        this.num_cols = data[0].length;
        this.shape = [this.num_rows, this.num_cols]
    }
    
    // okay in here we'll in theory handle two cases 
    // case 1: B is a scalar 
    // case 2: B is a vector/matrix where we'll broadcast the vector accross the matrix
    static elementWiseOperation(A, B, operation) {
        // Case 1: B is a scalar (e.g., matrix + scalar)
        if (typeof B === "number") {
            return new Matrix(A.data.map(row => row.map(value => operation(value, B))));
        } else if (B instanceof Matrix) {
            // Case 2: B is a matrix (check compatibility)
            // Case 2.1: A and B have the same number of rows, but different number of columns
            if (A.num_rows === B.num_rows && A.num_cols !== B.num_cols) {
                // broadcasting: assume B is a vector, broadcast it across the matrix
                // only going going to handle broadcasting alog the cols... so only care if num.cols == 1
                if (B.num_cols === 1) {
                    // okay so i dont need to use the i index for B here becuase it is only a vector
                    // like for example
                    // first we grab the all the dye emissions at the nth wavelength
                    // then we grab the ith dye at that nth_wavelength and use the nth_wavelength of the oprand dye 
                    // to do the arithmetic
                    return new Matrix(A.data.map((nth_row, n) => nth_row.map((ith_col, i) => operation(ith_col, B.data[n][0]))));
                } else {
                    // ill just make sure to use only compatible matricies along this dimension
                    throw new Error("Matrices are incompatible for broadcasting: number of columns in B must be 1 for broadcasting.");
                }
            }
            // Case 2.2: A and B have the same dimensions (element-wise operation)
            if (A.num_rows === B.num_rows && A.num_cols === B.num_cols) {
                
                // I'm kinda getting the hang of using the map and arrow functions... here is the same as the broadcasting case
                // except this time we are grabbing every element of B
                return new Matrix(A.data.map((ith_row, i) => ith_row.map((jth_col, j) => operation(jth_col, B.data[i][j]))));
            }
            // Case 2.3: A and B have different rows and columns (incompatible for element-wise operation)
            throw new Error("Matrices must have the same dimensions for element-wise operations.");
        } else {
            // just incase something else gets in there
            throw new Error("Second operand must be a scalar or a Matrix.");
        }
    }
    
    // for the operations we can use a simple arrow function
    add(B){
        return Matrix.elementWiseOperation(this, B, (a, b) => a + b)
    }

    subtract(B){
        return Matrix.elementWiseOperation(this, B, (a, b) => a - b)
    }

    multiply(B){
        return Matrix.elementWiseOperation(this, B, (a, b) => a * b)
    }

    divide(B){
        return Matrix.elementWiseOperation(this, B, (a, b) => a / b)
    }

    // basically A_ij -> A_ji
    transpose(){
        // gonna use for loops here intead of the map or array methods
        const A_transposed = [];
        // iterate through the columns
        for (let j = 0; j < this.num_cols; j++){
            const row = []; // make a new row for each column
            // iterate through the rows
            for (let i = 0; i < this.num_rows; i++){
                // this way for a single row i, it will iterate over j cols and add them to the new row array
                row.push(this.data[i][j]); 
            }
            // then send that row to the new array
            A_transposed.push(row)
        }
        return new Matrix(A_transposed)
    }

    dot(B){
        // TODO: Can make this cleaner 
        // I'm not sure this is the best way to do this when it comes to numerical stability
        // https://en.wikipedia.org/wiki/Dot_product#Algorithms seems like maybe not

        // check to make sure that B is an instance of a Matrix
        if (B instanceof Matrix){
            // need to check dimensions here first 
            if (this.num_cols != B.num_rows){
                throw new Error("Matrices are incompatible for dot product: A.num_cols must equal B.num_rows.");
            } else {
                // need a vector to store each line of the dot product thats computed
                const dot_result = [];
                for (let ith_row = 0; ith_row < this.num_rows; ith_row++){
                    // each row of the output matrix  corresponds to a row from the A matrix 
                    // so make a new vector here
                    let dot_row = [];
                    for (let jth_col = 0; jth_col < B.num_cols; jth_col++){
                        // now grab the column from B and iterate over the elements of both A_row and B_col
                        // each should have the same number of elements by definition 
                        let k_sum = 0;
                        for (let kth_element = 0; kth_element < this.num_cols; kth_element++){
                            // compute the sum
                            k_sum += this.data[ith_row][kth_element] * B.data[kth_element][jth_col];
                        }
                        // add the sum for that elements of that row-col pair as an entry in the row
                        dot_row.push(k_sum)
                    }
                    // then add the row as a output of the dot operation
                    dot_result.push(dot_row);
                }
                return new Matrix(dot_result)
            }
        } else {
            throw new Error("Scond operand must be Matrix.")
        }
    }

    // TODO: Matrix inversion using Gaussian Elimination
    invert() {
    }      

    //TODO only handles 2D matrices
    find_max(dimension){
        // search the rows of the matrix
        if (dimension == 0){
            // iterate through each row
            let output = [];
            for (let nth_row = 0; nth_row < this.num_rows; nth_row++){
                row = this.data[nth_row];
                row_max_val = Math.max(row);
                row_max_ind = row.indexOf(row_max_val)
            }
        }if (dimension == 1) {
            
        } else {
            console.error('Find max function only supports 2D matrices. Please choose from dimenions [0,1].');
        }
    }

    // a utility function to print the matrices
    toString() {
        return this.data.map(row => row.join("\t")).join("\n");
    }
}


const componentDataset = new Dataset();
const testDataset = new Dataset();

document.getElementById("processButton").addEventListener("click", async () => {
    
    const componentFile = document.getElementById("component_data").files[0];
    const testFile = document.getElementById('test_data').files[0];

    try {
        // Read both files
        const componentContents = await readFile(componentFile);
        const testContents = await readFile(testFile);

        // Load data into datasets
        componentDataset.loadFromContents(componentContents);
        testDataset.loadFromContents(testContents);

        // next we need to get the common wavelength range we want to use
        const wavel_upper = Math.max(componentDataset.wavel_max, testDataset.wavel_max);
        const wavel_lower = Math.min(componentDataset.wavel_min, testDataset.wavel_min);
        
        // we want 1 nm wavelength steps 
        const common_wavels = [];
        for (let nth_sample = wavel_lower; nth_sample <= wavel_upper; nth_sample++){
            common_wavels.push(nth_sample);
        }

        // next would be to interpolate the wavelengths gauranteeing they match 
        let interpolated_components = new Interpolator(componentDataset, common_wavels).interpolateDataset();
        let interpolated_test_data = new Interpolator(testDataset, common_wavels).interpolateDataset();

        const A = new Matrix(interpolated_components.spectra);
        const B = new Matrix(interpolated_test_data.spectra);

        // we want to interate over the column vectors of B and solve Ax = b
        // in each loop we want to solve for the coefficient vector of the least squares hyperplane
        // to do this solve x = (A^t A)^-1 A^t b 
        // see 
        let AtDotA = A.transpose().dot(A);
        console.log(AtDotA.invert())

    } catch (error) {
        console.error("Error reading files:", error);
    }
});

