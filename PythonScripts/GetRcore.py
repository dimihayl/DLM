from ROOT import TFile, TGraphAsymmErrors

def evaluate_tgraphs(input_file, mt):
  """
  Opens a ROOT file, retrieves TGraphs, evaluates them at a mass, and formats the output.

  Args:
      input_file: Path to the ROOT file.
      currently I tend to use ./BigLoop_v7_all/Results_pp.root
      mt: Mass value (float) for evaluation.

  Returns:
      A string containing the formatted output with central value and errors.
  """

  # Open the ROOT file
  root_file = TFile.Open(input_file)
  if not root_file:
    raise FileNotFoundError(f"Could not open file: {input_file}")

  # Get the TGraphs
  try:
    gOtonScript_original = root_file.Get("gOtonScript_original")
    gOtonScript_original_up = root_file.Get("gOtonScript_original_up")
    gOtonScript_original_low = root_file.Get("gOtonScript_original_low")
    gOtonScript_remake = root_file.Get("gOtonScript_remake")
    gOtonScript_remake_up = root_file.Get("gOtonScript_remake_up")
    gOtonScript_remake_low = root_file.Get("gOtonScript_remake_low")
  except KeyError:
    raise ValueError(f"TGraphs not found in file: {input_file}")

  # Close the ROOT file
  root_file.Close()

  # Evaluate TGraphs at the closest index
  original_value = gOtonScript_original.Eval(mt)
  original_yerr_up = gOtonScript_original_up.Eval(mt)-original_value
  original_yerr_low = original_value-gOtonScript_original_low.Eval(mt)
  original_avg_val = 0.5*(2*original_value+original_yerr_up-original_yerr_low)
  original_avg_err = 0.5*(original_yerr_up+original_yerr_low)

  remake_value = gOtonScript_remake.Eval(mt)
  remake_yerr_up = gOtonScript_remake_up.Eval(mt)-remake_value
  remake_yerr_low = remake_value-gOtonScript_remake_low.Eval(mt)
  remake_avg_val = 0.5*(2*remake_value+remake_yerr_up-remake_yerr_low)
  remake_avg_err = 0.5*(remake_yerr_up+remake_yerr_low)

  # Format output string
  print(f"\n@ mT = {mt:.2f} GeV")
  print(f"OLD r_core = {original_value:.3f} +{original_yerr_up:.3f} -{original_yerr_low:.3f} fm (min = {original_value-original_yerr_low:.3f}, max = {original_value+original_yerr_up:.3f})")
  print(f"    symmetric result: {original_avg_val:.3f} +/- {original_avg_err:.3f} fm")
  print(f"NEW r_core = {remake_value:.3f} +{remake_yerr_up:.3f} -{remake_yerr_low:.3f} fm (min = {remake_value-remake_yerr_low:.3f}, max = {remake_value+remake_yerr_up:.3f})")
  print(f"    symmetric result: {remake_avg_val:.3f} +/- {remake_avg_err:.3f} fm")
  print(f"")

# Example usage (assuming you have argparse installed)
import argparse

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Evaluate TGraphs in ROOT file")
  parser.add_argument("InputFile", help="Path to the ROOT file")
  parser.add_argument("mT", type=float, help="Mass value (GeV) for evaluation")
  args = parser.parse_args()

  try:
    output = evaluate_tgraphs(args.InputFile, args.mT)
    #print(output)
  except (FileNotFoundError, ValueError) as e:
    print(f"Error: {e}")