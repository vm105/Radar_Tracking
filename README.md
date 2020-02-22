- 2d cfar was implemented as follows:
  - an array of zeros was created of the size of RDM
  - the index of cell under test was determined per iteration
  - the indices of the left and right halves of the training matrix were determined
  - the values over the matrix was summed and averaged
  - the cfar threshold was determined using an offset threshold
  - values lesser than the cfar threshold were supressed

- Selection of training and guard cells:
  - this was done through observation iteratively
  - 3 training cells were selected for each half of the rows of the training matrix
  - 3 training cells were selected for each half of the cols of the training matrix
  - 2 guard cells were selected for each half of the rows of the training matrix
  - 2 guard  cells were selected for each half of the cols of the training matrix
  - The offset was chosen to be 5 dB


- Steps taken to suppress edge cells:
  - this was done by using a new array initialized with zeros
