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
        this.interpoalted_spectra = []
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
                const output_vector = [this.input_wavels[ind],this.spectra[ind],this.input_wavels[ind+1],this.spectra[ind+1],x_q];

                return output_vector
            }
        } 
    }

    _interpolatePoints(output_vector){

        // just going to unpack the output vector real quick
        const xi = output_vector[0];
        const Yi = output_vector[1];
        const xi_1 = output_vector[2];
        const Yi_1 = output_vector[3];
        const xq = output_vector[4];
        
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
            this.interpoalted_spectra.push(Yq);
        }
        let interpolated_dataset = new Dataset()
        
        return interpolated_dataset.loadFromVariables(this.query_points, this.interpoalted_spectra, this.headers)

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

        // next would be to interpolate the wavelengths
        const interpolated_components = new Interpolator(componentDataset, common_wavels).interpolateDataset();
        const interpolated_test_data = new Interpolator(testDataset, common_wavels).interpolateDataset();

        console.log(interpolated_components)

    } catch (error) {
        console.error("Error reading files:", error);
    }
});

